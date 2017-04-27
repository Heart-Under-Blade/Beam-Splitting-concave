#include "Particle.h"

Particle::Particle() {}

void Particle::Init(int facetNum, const complex &refrIndex, double symAngle,
					double size)
{
	m_facetNum = facetNum;
	m_refractionIndex = refrIndex;
	m_symmetryAngle = symAngle;
	m_mainSize = size;
}

void Particle::Rotate(double beta, double gamma, double alpha)
{
	SetRotateMatrix(beta, gamma, alpha);

	for (int i = 0; i < m_facetNum; ++i)
	{
		for (int j = 0; j < facets[i].polygon.size; ++j)
		{
			RotatePoint(defaultState.facets[i].polygon.arr[j],
						facets[i].polygon.arr[j]);
		}
	}

	RotateNormals();

	for (int i = 0; i < m_facetNum; ++i)
	{
		RotatePoint(defaultState.facets[i].center, facets[i].center);
	}
}

const double &Particle::GetMainSize() const
{
	return m_mainSize;
}

const complex &Particle::GetRefractionIndex() const
{
	return m_refractionIndex;
}

bool Particle::IsUnshadowedExternal(int facetId) const
{
	for (int i = 0; i < m_unshadowedExternalFacets.size; ++i)
	{
		if (facetId == m_unshadowedExternalFacets.arr[i])
		{
			return true;
		}
	}

	return false;
}

//REF неправильно названо (переименовать)
bool Particle::IsShadowedInternal(int facetId) const
{
	for (int i = 0; i < m_shadowedInternalFacets.size; ++i)
	{
		if (facetId == m_shadowedInternalFacets.arr[i])
		{
			return true;
		}
	}

	return false;
}

void Particle::SetDefaultNormals()
{
	for (int i = 0; i < m_facetNum; ++i)
	{
		defaultState.facets[i].SetNormal();
	}
}

void Particle::SetActualState()
{
	for (int i = 0; i < m_facetNum; ++i)
	{
		facets[i] = defaultState.facets[i];
	}
}

const double &Particle::GetSymmetryAngle() const
{
    return m_symmetryAngle;
}

void Particle::SetRotateMatrix(double beta, double gamma, double alpha)
{
	double cosA, cosB, cosG,
			sinA, sinB, sinG;

	sincos(alpha, &sinA, &cosA);
	sincos(beta,  &sinB, &cosB);
	sincos(gamma, &sinG, &cosG);

	double cosAcosB = cosA*cosB;
	double sinAcosG = sinA*cosG;
	double sinAsinG = sinA*sinG;

	m_rotMatrix[0][0] = cosAcosB*cosG - sinAsinG;
	m_rotMatrix[1][0] = sinAcosG*cosB + cosA*sinG;
	m_rotMatrix[2][0] = -sinB*cosG;

	m_rotMatrix[0][1] = -(cosAcosB*sinG + sinAcosG);
	m_rotMatrix[1][1] = cosA*cosG - sinAsinG*cosB;
	m_rotMatrix[2][1] = sinB*sinG;

	m_rotMatrix[0][2] = cosA*sinB;
	m_rotMatrix[1][2] = sinA*sinB;
	m_rotMatrix[2][2] = cosB;
}

void Particle::RotateNormals()
{
	for (int i = 0; i < m_facetNum; ++i)
	{
		RotatePoint(defaultState.facets[i].in_normal, facets[i].in_normal);
	}

	SetDParams();

	for (int i = 0; i < m_facetNum; ++i)
	{
		facets[i].ex_normal = -facets[i].in_normal;
		facets[i].ex_normal.d_param = -facets[i].in_normal.d_param;
	}
}

void Particle::SetDParams()
{
	for (int i = 0; i < m_facetNum; ++i)
	{
		double d = DotProduct(facets[i].polygon.arr[0], facets[i].in_normal);
		facets[i].in_normal.d_param = -d;
	}
}

void Particle::RotatePoint(const Point3f &point, Point3f &result)
{
	result.cx = point.cx*m_rotMatrix[0][0] + point.cy*m_rotMatrix[0][1] + point.cz*m_rotMatrix[0][2];
	result.cy = point.cx*m_rotMatrix[1][0] + point.cy*m_rotMatrix[1][1] + point.cz*m_rotMatrix[1][2];
	result.cz = point.cx*m_rotMatrix[2][0] + point.cy*m_rotMatrix[2][1] + point.cz*m_rotMatrix[2][2];
}
