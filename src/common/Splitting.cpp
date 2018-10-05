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

void Splitting::ComputeSplittingParams(const Point3f &dir, const Point3f &normal,
									   bool isInside)
{
	m_normal = normal;
	r = dir/cosA - m_normal;

	if (isInside)
	{
		reRiEff = ComputeEffectiveReRi();
		s = 1.0/(reRiEff*cosA*cosA) - Point3f::Norm(r);
	}
	else
	{
		double cosA2 = cosA * cosA;
		double tmp = m_cRiRe + cosA2 - 1.0;

		if (m_cRiIm > FLT_EPSILON)
		{
			tmp = sqrt(tmp*tmp + m_cRiIm);
		}

		tmp = (m_cRiRe + 1.0 - cosA2 + tmp)/2.0;
		s = (tmp/cosA2) - Point3f::Norm(r);
	}
}

double Splitting::ComputeEffectiveReRi() const
{
	return (m_cRiRe + sqrt(m_cRiRe2 + m_cRiIm/(cosA*cosA)))/2.0;
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

void Splitting::SetBeams(const Polygon &beamShape)
{
	inBeam.Clear();
	inBeam.SetPolygon(beamShape);

	outBeam.Clear();
	outBeam.SetPolygon(beamShape);
}

void Splitting::SetNormal(const Point3f &normal)
{
	m_normal = normal;
}

double Splitting::ComputeSegmentOpticalPath(const Beam &beam, const Point3f &facetPoint) const
{
	double tmp = Point3d::DotProduct(beam.direction, facetPoint);
	double path = fabs(tmp + beam.front); // refractive index of external media = 1

#ifdef _DEBUG // DEB
//	if (tmp + beam.front < 0)
//		int fff = 0;
	Point3f dd = beam.Center();
	Point3f pd = beam.Center() + (beam.direction * 10);

	/* ПРоверка на нахождения плоскости за пучком
	 * REF: вынести в случай невыпуклых частиц, т.к. характерно только для них */
	double cosB = Point3d::DotProduct(Point3d::Normalize(facetPoint - beam.Center()), beam.direction);
#else
	Vector3f vectorToCenter = facetPoint - beam.Center();
	Point3f::Normalize(vectorToCenter);
	double cosB = Point3f::DotProduct(vectorToCenter, beam.direction);
#endif

	if (cosB < 0 && beam.isInside)
	{
		path = -path;
	}

	if (beam.isInside)
	{
		path *= sqrt(reRiEff);
	}

	return path;
}

void Splitting::ComputeCosA(const Point3f &normal, const Point3f &incidentDir)
{
	cosA = Point3f::DotProduct(normal, incidentDir);
}

double Splitting::ComputeIncidentOpticalPath(const Point3f &direction,
											 const Point3f &facetPoint)
{
	return FAR_ZONE_DISTANCE + Point3f::DotProduct(direction, facetPoint);
}

double Splitting::ComputeOutgoingOpticalPath(const Beam &beam)
{
	return FAR_ZONE_DISTANCE + beam.front;
}

complex Splitting::GetRi() const
{
	return m_ri;
}
