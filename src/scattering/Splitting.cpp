#include "Splitting.h"

#include <float.h>

#define EPS_ZERO	1.7453292519943295769148298069306e-10	//cos(89.99999999)
#define EPS_ONE		0.99999999998254670756866631966593		//1 - cos(89.99999999)

Splitting::Splitting(const complex &ri)
{
	m_ri = ri;
	double re = real(m_ri);
	double im = imag(m_ri);
	m_cRi0 = re*re - im*im;
	m_cRi1 = m_cRi0 * m_cRi0;
	m_cRi2 = 4*re*re*im;
}

void Splitting::ComputeSplittingParams(const Vector3f &dir, const Vector3f &normal,
									   bool isInside)
{
	m_hasOutBeam = true;
	cosA = Point3f::DotProduct(normal, dir);

	if (cosA > EPS_ONE)
	{
		facetNormal = normal;
		type = IncidenceType::Normal;
	}
	else
	{
		cosA2 = cosA * cosA;
		r = dir/cosA - normal;

		if (isInside)
		{
			reRiEff = ComputeEffectiveReRi(cosA2);
			s = 1.0/(reRiEff*cosA2) - Point3f::Norm(r);

			if (s < DBL_EPSILON)
			{
				m_hasOutBeam = false;
				type = IncidenceType::CompleteReflection;
			}
			else
			{
				facetNormal = normal;
				type = IncidenceType::Regular;
			}
		}
		else
		{
			facetNormal = normal;
			type = IncidenceType::Regular;

			double tmp = m_cRi0 + cosA2 - 1.0;

			if (m_cRi2 > FLT_EPSILON)
			{
				tmp = sqrt(tmp*tmp + m_cRi2);
			}

			tmp = (m_cRi0 + 1.0 - cosA2 + tmp)/2.0;
			s = (tmp/cosA2) - Point3f::Norm(r);
		}
	}
}

IncidenceType Splitting::GetIncidenceType() const
{
	return type;
}

double Splitting::ComputeEffectiveReRi(double cosA2) const
{
	return (m_cRi0 + sqrt(m_cRi1 + m_cRi2/cosA2))/2.0;
}

void Splitting::SetBeams(const Polygon &beamShape)
{
	beams.internal.Clear();
	beams.internal.SetPolygon(beamShape);

	beams.external.Clear();
	beams.external.SetPolygon(beamShape);

#ifdef _DEBUG // DEB
	beams.internal.pols.push_back(beamShape);
	beams.external.pols.push_back(beamShape);
#endif
}

bool Splitting::HasOutBeam()
{
	return m_hasOutBeam;
}

complex Splitting::GetRi() const
{
	return m_ri;
}
