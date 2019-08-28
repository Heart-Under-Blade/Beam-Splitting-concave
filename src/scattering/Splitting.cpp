#include "Splitting.h"

#include <float.h>

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
	m_r = dir/cosA - normal;
	reRiEff = ComputeEffectiveReRi();
	s = 1.0/(reRiEff*cosA*cosA) - Norm(m_r);
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
	double tmp = DotProductD(beam.direction, facetPoint);

	double path = fabs(tmp + beam.front); // refractive index of external media = 1

	/* ПРоверка на нахождения плоскости за пучком
	 * REF: вынести в случай невыпуклых частиц, т.к. характерно только для них */
	double cosB = DotProductD(NormalizeD(facetPoint - beam.Center()), beam.direction);

	if (cosB < 0)
	{
		path = -path;
	}

	if (beam.location == Location::In)
	{
		path *= sqrt(reRiEff);
	}

	return path;
}

// REF: создать отдельные классы RegularSplitting, CRSplitting, NormalSplitting

void Splitting::ComputeCRBeamParams(const Point3f &normal, const Beam &incidentBeam,
									Beam &inBeam)
{
	Point3f reflDir = m_r - normal;
	Normalize(reflDir);
	inBeam.SetLight(reflDir, incidentBeam.polarizationBasis);

	complex cv, ch;
	ComputeCRJonesParams(cv, ch);

	inBeam.J = incidentBeam.J;
	inBeam.MultiplyJonesMatrix(cv, ch);

	if (m_isOpticalPath)
	{
		double path = ComputeSegmentOpticalPath(incidentBeam, inBeam.Center());
#ifndef _DEBUG // DEB
//		inBeam.ops = incidentBeam.ops;
//		inBeam.ops.push_back(path);
#endif
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

void Splitting::RefractIn(const Vector3f &r, const Vector3f &normal,
						  Vector3f &newDir)
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
	newDir = (r/tmp) - normal;
	Normalize(newDir);
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
	Point3f reflDir = m_r - normal;
	Normalize(reflDir);
	inBeam.SetLight(reflDir, incidentBeam.polarizationBasis);

	Point3f refrDir = m_r/sqrt(s) + normal;
	Normalize(refrDir);
	outBeam.SetLight(refrDir, incidentBeam.polarizationBasis);

	ComputeRegularJonesParams(normal, incidentBeam, inBeam, outBeam);

	if (m_isOpticalPath)
	{
		double path = ComputeSegmentOpticalPath(incidentBeam, inBeam.Center());
#ifndef _DEBUG // DEB
//		inBeam.ops = incidentBeam.ops;
//		outBeam.ops = incidentBeam.ops;
//		outBeam.ops.push_back(path);
#endif
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

//	if (m_isOpticalPath)
	{
		double path = ComputeSegmentOpticalPath(incidentBeam, inBeam.Center());
		path += incidentBeam.opticalPath;
#ifndef _DEBUG // DEB
//		inBeam.ops = incidentBeam.ops;
//		inBeam.ops.push_back(path);
//		outBeam.ops = incidentBeam.ops;
//		outBeam.ops.push_back(path);
#endif
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

	RefractIn(r, -facetNormal, inBeam.direction);

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

void Splitting::ReflectExternal(const Vector3f &oldDir, const Vector3f &normal,
						   Vector3f &newDir)
{
	Vector3f r = normal + oldDir/cosA;
	newDir = normal + r;
	Normalize(newDir);
}

void Splitting::ReflectInternal(const Vector3f &oldDir, const Vector3f &normal,
								Vector3f &newDir)
{
	ComputeCosA(oldDir, normal);
	ComputeSplittingParams(oldDir, normal);
	newDir = m_r - normal;
	Normalize(newDir);
}

void Splitting::RefractOut(const Vector3f &oldDir, const Vector3f &normal, Vector3f &newDir)
{
	ComputeCosA(oldDir, normal);
	ComputeSplittingParams(oldDir, normal);
	newDir = m_r/sqrt(s) + normal;
	Normalize(newDir);
}

Point3f Splitting::ChangeBeamDirectionConvex(const Vector3f &oldDir,
											 const Vector3f &normal,
											 Location loc)
{
	Point3f newDir;

	if (loc == Location::Out)
	{
		ComputeCosA(oldDir, -normal);
		Vector3f r = normal + oldDir/cosA;
		RefractIn(r, normal, newDir);
	}
	else
	{
		ReflectInternal(oldDir, normal, newDir);
	}

	return newDir;
}

Point3f Splitting::ChangeBeamDirection(const Vector3f &oldDir,
									   const Vector3f &normal,
									   Location oldLoc, Location loc)
{
	Point3f newDir;

	if (oldLoc == loc)
	{
		if (oldLoc == Location::Out)
		{
			ComputeCosA(oldDir, -normal);
			ReflectExternal(oldDir, normal, newDir);
		}
		else
		{
			ReflectInternal(oldDir, normal, newDir);
		}
	}
	else
	{
		if (oldLoc == Location::Out)
		{
			ComputeCosA(oldDir, -normal);
			Vector3f r = normal + oldDir/cosA;
			RefractIn(r, normal, newDir);
		}
		else
		{
			RefractOut(oldDir, normal, newDir);
		}
	}

	return newDir;
}

complex Splitting::GetRi() const
{
	return m_ri;
}

double Splitting::ComputeEffectiveReRi() const
{
	return (m_cRiRe + sqrt(m_cRiRe2 + m_cRiIm/(cosA*cosA)))/2.0;
}
