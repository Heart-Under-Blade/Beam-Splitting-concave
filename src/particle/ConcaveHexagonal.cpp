#include "ConcaveHexagonal.h"

ConcaveHexagonal::ConcaveHexagonal(double cavityDept)
{
	m_cavityDept = cavityDept;
}

ConcaveHexagonal::ConcaveHexagonal(double radius, double halfHeight, const complex &refractionIndex,
								   double cavityDept)
	: Hexagonal(radius, halfHeight, refractionIndex)
{
	m_cavityDept = cavityDept;

	SetFacetParams();
	SetBaseFacets();
	SetOriginCavityPoints();

	SetCavities(m_originBases.top, m_originBases.bottom, m_defaultCavities);
	SetSideFacets(m_originBases.top, m_originBases.bottom,
				  BASE_VERTEX_NUM, 2*BASE_VERTEX_NUM);

	SetDefaultNormals();
	SetActualNormals();
	SetCenters();
}

void ConcaveHexagonal::Rotate(double beta, double gamma, double alpha)
{
	Particle::Rotate(beta, gamma, alpha);

	Point3f baseTop[BASE_VERTEX_NUM], baseBottom[BASE_VERTEX_NUM];
	RotateBaseFacets(baseTop, baseBottom);

	SetSideFacets(baseTop, baseBottom, BASE_VERTEX_NUM, 2*BASE_VERTEX_NUM);

	CavityPoints cavities;
	RotatePoint(m_defaultCavities.top, cavities.top);
	RotatePoint(m_defaultCavities.bottom, cavities.bottom);
	SetCavities(baseTop, baseBottom, cavities);

	RotateNormals();

	for (int i = 0; i < facetNum; ++i)
	{
		RotatePoint(m_originCenters[i], centers[i]);
	}
}

void ConcaveHexagonal::SetCenters()
{
	for (int i = 0; i < facetNum; ++i)
	{
		m_originCenters[i] = facets[i].polygon.Center();
		centers[i] = m_originCenters[i];
	}
}

void ConcaveHexagonal::SetFacetParams()
{
	facetNum = BASE_VERTEX_NUM + BASE_FACET_NUM * BASE_VERTEX_NUM;

	// top facet (triangles)
	for (int i = 0; i < BASE_VERTEX_NUM; ++i)
	{
		facets[i].polygon.size = CAVITY_FACET_VERTEX_NUM;
	}

	// side facet
	for (int i = BASE_VERTEX_NUM; i < 2*BASE_VERTEX_NUM; ++i)
	{
		facets[i].polygon.size = SIDE_VERTEX_NUMBER;
	}

	// bottom facet (triangles)
	for (int i = 2*BASE_VERTEX_NUM; i < 3*BASE_VERTEX_NUM; ++i)
	{
		facets[i].polygon.size = CAVITY_FACET_VERTEX_NUM;
	}

	// shodowed and unshadowed
	for (int i = BASE_VERTEX_NUM; i < 2*BASE_VERTEX_NUM; ++i)
	{
		m_unshadowedExternalFacets.Add(i);
		m_shadowedInternalFacets.Add(i);
	}
}

void ConcaveHexagonal::SetCavities(const Point3f *baseTop, const Point3f *baseBottom,
								   const CavityPoints &cavities)
{
	SetCavityFacets(0, BASE_VERTEX_NUM, baseTop, cavities.top); // top facets (triangles)
	SetCavityFacets(2*BASE_VERTEX_NUM, 3*BASE_VERTEX_NUM, baseBottom, cavities.bottom); // bottom facets (triangles)
}

void ConcaveHexagonal::SetCavityFacets(int start, int end,
									   const Point3f *baseFacet,
									   const Point3f &cavityPoint)
{
	// base facet point indices
	int p0 = BASE_VERTEX_NUM-1;
	int p1 = 0;

	for (int i = start; i < end; ++i)
	{
		facets[i].polygon.arr[0] = baseFacet[p0];
		facets[i].polygon.arr[1] = baseFacet[p1];
		facets[i].polygon.arr[2] = cavityPoint;
		p0 = p1;
		++p1;
	}
}

void ConcaveHexagonal::SetOriginCavityPoints()
{
	m_defaultCavities.top = Point3f(0, 0, m_halfHeight - m_cavityDept);
	m_defaultCavities.bottom = Point3f(0, 0, -m_halfHeight + m_cavityDept);
}

void ConcaveHexagonal::SetBaseNormals()
{
	// top cavity facets
	for (int i = 0; i < BASE_VERTEX_NUM; ++i)
	{
		m_originNormals[i] = -facets[i].polygon.Normal();
	}

	// bottom cavity facets
	for (int i = 2*BASE_VERTEX_NUM; i < 3*BASE_VERTEX_NUM; ++i)
	{
		m_originNormals[i] = -facets[i].polygon.Normal();
	}
}

void ConcaveHexagonal::SetDefaultNormals()
{
	SetBaseNormals();
	SetSideNormals(6);
}
