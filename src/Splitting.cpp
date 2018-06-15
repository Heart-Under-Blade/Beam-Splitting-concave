#include "Splitting.h"

#include <float.h>

#define EPS_COS_90		1.7453292519943295769148298069306e-10	//cos(89.99999999)
#define EPS_COS_00		0.99999999998254670756866631966593		//1 - cos(89.99999999)

Splitting::Splitting(bool isOpticalPath)
{
	m_isOpticalPath = isOpticalPath;
}

void Splitting::ComputeRiParams(const complex &ri)
{
	m_ri = ri;
	double re = real(m_ri);
	double im = imag(m_ri);
	m_cRiRe = re*re - im*im;
	m_cRiRe2 = m_cRiRe * m_cRiRe;
	m_cRiIm = 4*re*re*im;
}

void Splitting::ComputeSplittingParams(const Point3f &dir, const Point3f &normal)
{
	r = dir/cosA - normal;
	reRiEff = ComputeEffectiveReRi();
	s = 1.0/(reRiEff*cosA*cosA) - Norm(r);
}

bool Splitting::IsCompleteReflection()
{
	return s <= DBL_EPSILON;
}

bool Splitting::IsNormalIncidence()
{
	return cosA >= EPS_COS_00;
}

bool Splitting::IsIncident()
{
	return cosA >= EPS_COS_90;
}

double Splitting::ComputeSegmentOpticalPath(const Beam &beam, const Point3f &facetPoint) const
{
	double tmp = DotProduct(beam.direction, facetPoint);
	double path = fabs(tmp + beam.front); // refractive index of external environment = 1

	if (beam.location == Location::In)
	{
		path *= sqrt(reRiEff);
	}

	return path;
}

void Splitting::ComputeCRBeamParams(const Point3f &normal, const Beam &incidentBeam,
									Beam &inBeam)
{
	Point3f reflDir = r - normal;
	Normalize(reflDir);
	inBeam.SetLight(reflDir, incidentBeam.polarizationBasis);

	complex cv, ch;
	ComputeCRJonesParams(cv, ch);

	inBeam.J = incidentBeam.J;
	inBeam.MultiplyJonesMatrix(cv, ch);

	if (m_isOpticalPath)
	{
		double path = ComputeSegmentOpticalPath(incidentBeam, inBeam.Center());
		path += incidentBeam.opticalPath;
		inBeam.AddOpticalPath(path);
	}
}

void Splitting::ComputeRegularJonesParams(const Point3f &normal, const Beam &incidentBeam, Beam &inBeam, Beam &outBeam)
{
	inBeam.J = incidentBeam.J;
	outBeam.J = incidentBeam.J;

	double cosG = DotProduct(normal, outBeam.direction);

	complex tmp0 = m_ri * cosA;
	complex tmp1 = m_ri * cosG;
	complex Tv0 = tmp1 + cosA;
	complex Th0 = tmp0 + cosG;

	complex tmp = 2.0 * tmp0;
	outBeam.MultiplyJonesMatrix(tmp/Tv0, tmp/Th0);

	complex Tv = cosA - tmp1;
	complex Th = tmp0 - cosG;
	inBeam.MultiplyJonesMatrix(Tv/Tv0, Th/Th0);
}

void Splitting::ComputeInternalRefractiveDirection(const Vector3f &r,
												   const Vector3f &normal,
												   Vector3f &dir)
{
	double cosA2 = cosA * cosA;
	double tmp = m_cRiRe + cosA2 - 1.0;

	if (m_cRiIm > FLT_EPSILON)
	{
		tmp = sqrt(tmp*tmp + m_cRiIm);
	}

	tmp = (m_cRiRe + 1.0 - cosA2 + tmp)/2.0;
	tmp = (tmp/cosA2);
	tmp -= Norm(r);
	tmp = sqrt(tmp);
	dir = (r/tmp) - normal;
	Normalize(dir);
}

void Splitting::ComputeCRJonesParams(complex &cv, complex &ch)
{
	const double bf = reRiEff*(1.0 - cosA*cosA) - 1.0;
	double im = (bf > 0) ? sqrt(bf) : 0;

	const complex sq(0, im);
	complex tmp0 = m_ri * cosA;
	complex tmp1 = m_ri * sq;

	cv = (cosA - tmp1)/(tmp1 + cosA);
	ch = (tmp0 - sq)/(tmp0 + sq);
}

void Splitting::ComputeCosA(const Point3f &normal, const Point3f &incidentDir)
{
	cosA = DotProduct(normal, incidentDir);
}

void Splitting::ComputeRegularBeamsParams(const Point3f &normal,
										  const Beam &incidentBeam,
										  Beam &inBeam, Beam &outBeam)
{
	Point3f reflDir = r - normal;
	Normalize(reflDir);
	inBeam.SetLight(reflDir, incidentBeam.polarizationBasis);

	Point3f refrDir = r/sqrt(s) + normal;
	Normalize(refrDir);
	outBeam.SetLight(refrDir, incidentBeam.polarizationBasis);

	ComputeRegularJonesParams(normal, incidentBeam, inBeam, outBeam);

	if (m_isOpticalPath)
	{
		double path = ComputeSegmentOpticalPath(incidentBeam, inBeam.Center());
		path += incidentBeam.opticalPath;
		inBeam.AddOpticalPath(path);
		outBeam.AddOpticalPath(path);
	}
}

void Splitting::ComputeNormalBeamParams(const Beam &incidentBeam,
										Beam &inBeam, Beam &outBeam)
{
	const Point3f &dir = incidentBeam.direction;
	inBeam.SetLight(-dir, incidentBeam.polarizationBasis);
	outBeam.SetLight(dir, incidentBeam.polarizationBasis);

	inBeam.J = incidentBeam.J;
	outBeam.J = incidentBeam.J;

	complex temp;

	temp = (2.0 * m_ri)/(1.0 + m_ri); // OPT: вынести целиком
	outBeam.MultiplyJonesMatrix(temp, temp);

	temp = (1.0 - m_ri)/(1.0 + m_ri); // OPT: вынести целиком
	inBeam.MultiplyJonesMatrix(temp, -temp);

	if (m_isOpticalPath)
	{
		double path = ComputeSegmentOpticalPath(incidentBeam, inBeam.Center());
		path += incidentBeam.opticalPath;
		inBeam.AddOpticalPath(path);
		outBeam.AddOpticalPath(path);
	}
}

void Splitting::ComputeNormalBeamParamsExternal(const Light &incidentLight,
												Beam &inBeam, Beam &outBeam)
{
	inBeam.SetLight(incidentLight);
	outBeam.SetLight(-incidentLight.direction, incidentLight.polarizationBasis);

	inBeam.J.m11 = 2.0/(m_ri + 1.0); // OPT: вынести
	inBeam.J.m22 = inBeam.J.m11;

	outBeam.J.m11 = (m_ri - 1.0)/(m_ri + 1.0);
	outBeam.J.m22 = -outBeam.J.m11;
}

void Splitting::ComputeRegularBeamParamsExternal(const Point3f &facetNormal,
												 Beam &incidentBeam,
												 Beam &inBeam, Beam &outBeam)
{
	inBeam.polarizationBasis = incidentBeam.polarizationBasis;

	inBeam.J = incidentBeam.J;
	outBeam.J = incidentBeam.J;

	Point3f refrDir;
	const Point3f &dir = incidentBeam.direction;

	Point3f r = dir/cosA - facetNormal;

	refrDir = r - facetNormal;
	Normalize(refrDir);

	ComputeInternalRefractiveDirection(r, -facetNormal, inBeam.direction);

	outBeam.SetLight(refrDir, incidentBeam.polarizationBasis);

	double cosB = DotProduct(facetNormal, inBeam.direction);

	complex tmp0 = m_ri * cosA;
	complex tmp1 = m_ri * cosB;

	complex Tv0 = tmp0 + cosB;
	complex Th0 = tmp1 + cosA;

	complex Tv = tmp0 - cosB;
	complex Th = cosA - tmp1;
	outBeam.MultiplyJonesMatrix(Tv/Tv0, Th/Th0);

	double cos2A = 2.0*cosA;
	inBeam.MultiplyJonesMatrix(cos2A/Tv0, cos2A/Th0);
}

double Splitting::ComputeIncidentOpticalPath(const Point3f &direction,
											 const Point3f &facetPoint)
{
	return FAR_ZONE_DISTANCE + DotProduct(direction, facetPoint);
}

double Splitting::ComputeOutgoingOpticalPath(const Beam &beam)
{
	return FAR_ZONE_DISTANCE + beam.front;
}

Point3f Splitting::ChangeBeamDirection(const Vector3f &oldDir,
									   const Vector3f &normal,
									   Location oldLoc, Location loc)
{
	Point3f newDir;

	if (oldLoc == Location::Out) // refraction
	{
		ComputeCosA(oldDir, -normal);
		Vector3f r = normal + oldDir/cosA;

		if (oldLoc == loc)
		{
			newDir = normal + r;
			Normalize(newDir);
		}
		else
		{
			ComputeInternalRefractiveDirection(r, normal, newDir);
		}
	}
	else // reflection
	{
		ComputeCosA(oldDir, normal);
		ComputeSplittingParams(oldDir, normal);

		if (oldLoc == loc)
		{
			newDir = r - normal;
		}
		else
		{
			newDir = r/sqrt(s) + normal;
		}

		Normalize(newDir);
	}

	return newDir;
}

double Splitting::ComputeEffectiveReRi() const
{
	return (m_cRiRe + sqrt(m_cRiRe2 + m_cRiIm/(cosA*cosA)))/2.0;
}
