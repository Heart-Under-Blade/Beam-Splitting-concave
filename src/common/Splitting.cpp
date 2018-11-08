#include "Splitting.h"

#include <float.h>

#define EPS_COS_90		1.7453292519943295769148298069306e-10	//cos(89.99999999)
#define EPS_COS_00		0.99999999998254670756866631966593		//1 - cos(89.99999999)

Splitting::Splitting(const complex &ri)
{
	m_ri = ri;
	double re = real(m_ri);
	double im = imag(m_ri);
	m_cRiRe = re*re - im*im;
	m_cRiRe2 = m_cRiRe * m_cRiRe;
	m_cRiIm = 4*re*re*im;
}

void Splitting::ComputeParams(const Vector3f &dir, const Vector3f &normal,
							  bool isInside)
{
	m_hasOutBeam = true;
	m_normal = normal;
	cosA = Point3f::DotProduct(m_normal, dir);

	if (cosA > EPS_COS_00)
	{
		m_incidence = &m_normalIncidence;
	}
	else
	{
		r = dir/cosA - m_normal;

		if (isInside)
		{
			reRiEff = ComputeEffectiveReRi();
			s = 1.0/(reRiEff*cosA*cosA) - Point3f::Norm(r);

			if (s < DBL_EPSILON)
			{
				m_hasOutBeam = false;
				m_incidence = &m_completeReflectionIncidence;
			}
			else
			{
				m_incidence = &m_regularIncidence;
			}
		}
		else
		{
			m_incidence = &m_regularIncidence;

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
}

double Splitting::ComputeEffectiveReRi() const
{
	return (m_cRiRe + sqrt(m_cRiRe2 + m_cRiIm/(cosA*cosA)))/2.0;
}

void Splitting::SetNormal(const Point3f &normal)
{
	m_normal = normal;
}

bool Splitting::HasOutBeam()
{
	return m_hasOutBeam;
}

complex Splitting::GetRi() const
{
	return m_ri;
}

Incidence *Splitting::GetIncidence() const
{
	return m_incidence;
}
