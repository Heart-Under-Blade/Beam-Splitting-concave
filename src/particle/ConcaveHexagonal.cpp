#include "ConcaveHexagonal.h"
#include "global.h"
#include <algorithm>

ConcaveHexagonal::ConcaveHexagonal(const complex &refrIndex,
								   double diameter, double height,
								   double cavityDept)
{
	SetSize(diameter, height);
	m_cavityDept = cavityDept;

	double size = std::max(m_height, m_diameter);
	Init(18, refrIndex, M_PI/3, size);

	SetFacetParams();

	Facet baseTop, baseBottom;
	SetBases(baseTop, baseBottom);

	SetOriginCavityPoints();

	SetCavities(baseTop, baseBottom, m_defaultStateCavities);
	SetSides(baseTop, baseBottom);

	SetDefaultNormals();
	SetActualState();
	SetCenters();
}

void ConcaveHexagonal::SetCenters()
{
	for (int i = 0; i < m_facetNum; ++i)
	{
		defaultState.facets[i].SetCenter();
	}
}

void ConcaveHexagonal::SetFacetParams()
{
	SetSideFacetParams(BASE_VERTEX_NUM, 2*BASE_VERTEX_NUM);

	// top facet (triangles)
	for (int i = 0; i < BASE_VERTEX_NUM; ++i)
	{
		defaultState.facets[i].polygon.size = CAVITY_FACET_VERTEX_NUM;
	}

	// bottom facet (triangles)
	for (int i = 2*BASE_VERTEX_NUM; i < 3*BASE_VERTEX_NUM; ++i)
	{
		defaultState.facets[i].polygon.size = CAVITY_FACET_VERTEX_NUM;
	}

	// shodowed and unshadowed
	for (int i = BASE_VERTEX_NUM; i < 2*BASE_VERTEX_NUM; ++i)
	{
		m_unshadowedExternalFacets.Add(i);
		m_shadowedInternalFacets.Add(i);
	}
}

void ConcaveHexagonal::SetCavities(Facet &baseTop, Facet &baseBottom,
								   const CavityPoints &cavities)
{
	Point3f *bTop = baseTop.polygon.arr;
	Point3f *bBot = baseBottom.polygon.arr;

	SetCavityFacets(0, BASE_VERTEX_NUM, bTop, cavities.top); // top facets (triangles)
	SetCavityFacets(2*BASE_VERTEX_NUM, 3*BASE_VERTEX_NUM, bBot, cavities.bottom); // bottom facets (triangles)
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
		defaultState.facets[i].polygon.arr[0] = baseFacet[p0];
		defaultState.facets[i].polygon.arr[1] = baseFacet[p1];
		defaultState.facets[i].polygon.arr[2] = cavityPoint;
		p0 = p1;
		++p1;
	}
}

void ConcaveHexagonal::SetOriginCavityPoints()
{
	m_defaultStateCavities.top	  = Point3f(0, 0,  m_height/2 - m_cavityDept);
	m_defaultStateCavities.bottom = Point3f(0, 0, -m_height/2 + m_cavityDept);
}
