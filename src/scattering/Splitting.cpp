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
	cosA = Point3f::DotProduct(normal, dir);
	facetNormal = normal;

	if (cosA > EPS_ONE)
	{
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
				type = IncidenceType::TotalReflection;
			}
			else
			{
				type = IncidenceType::Regular;
			}
		}
		else
		{
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

void Splitting::ComputeSplittingParams(const Point3f &dir, const Point3f &normal)
{
	r = dir/cosA - normal;
	reRiEff = ComputeEffectiveReRi(cosA);
	s = 1.0/(reRiEff*cosA*cosA) - Point3f::Norm(r);
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
	beams.internal.SetDefault();
	beams.internal.SetPolygon(beamShape);

	beams.external.SetDefault();
	beams.external.SetPolygon(beamShape);
}

complex Splitting::GetRi() const
{
	return m_ri;
}

Point3f Splitting::ChangeBeamDirection(const Vector3f &oldDir,
									   const Vector3f &normal,
									   bool isInOld, bool isInNew)
{
	Point3f newDir;

	if (isInOld == isInNew)
	{
		if (!isInOld)
		{
			cosA = Point3f::DotProduct(oldDir, -normal);
			ReflectExternal(oldDir, normal, newDir);
		}
		else
		{
			ReflectInternal(oldDir, normal, newDir);
		}
	}
	else
	{
		if (!isInOld)
		{
			cosA = Point3f::DotProduct(oldDir, -normal);
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

Point3f Splitting::ChangeBeamDirectionConvex(const Vector3f &oldDir,
											 const Vector3f &normal, bool isIn)
{
	Point3f newDir;

	if (!isIn)
	{
		cosA = Point3f::DotProduct(oldDir, -normal);
		Vector3f r = normal + oldDir/cosA;
		RefractIn(r, normal, newDir);
	}
	else
	{
		ReflectInternal(oldDir, normal, newDir);
	}

	return newDir;
}

void Splitting::RefractIn(const Vector3f &r, const Vector3f &normal, Vector3f &newDir)
{
	double cosA2 = cosA * cosA;
	double tmp = m_cRi0 + cosA2 - 1.0;

	if (m_cRi2 > FLT_EPSILON)
	{
		tmp = sqrt(tmp*tmp + m_cRi2);
	}

	tmp = (m_cRi0 + 1.0 - cosA2 + tmp)/2.0;
	tmp = (tmp/cosA2);
	tmp -= Point3f::Norm(r);
	tmp = sqrt(tmp);
	newDir = (r/tmp) - normal;
	Point3f::Normalize(newDir);
}

void Splitting::ReflectExternal(const Vector3f &oldDir, const Vector3f &normal,
								Vector3f &newDir)
{
	Vector3f r = normal + oldDir/cosA;
	newDir = normal + r;
	Point3f::Normalize(newDir);
}

void Splitting::ReflectInternal(const Vector3f &oldDir, const Vector3f &normal, Vector3f &newDir)
{
	cosA = Point3f::DotProduct(oldDir, normal);
	ComputeSplittingParams(oldDir, normal);
	newDir = r - normal;
	Point3f::Normalize(newDir);
}

void Splitting::RefractOut(const Vector3f &oldDir, const Vector3f &normal, Vector3f &newDir)
{
	cosA = Point3f::DotProduct(oldDir, normal);
	ComputeSplittingParams(oldDir, normal);
	newDir = r/sqrt(s) + normal;
	Point3f::Normalize(newDir);
}
