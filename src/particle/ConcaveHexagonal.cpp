#include "ConcaveHexagonal.h"
#include "global.h"

ConcaveHexagonal::ConcaveHexagonal(double radius, double halfHeight, const complex &refractionIndex,
								   double cavityDept)
{
	m_symmetryAngle = M_PI/3;
	Particle::Init(radius, halfHeight, refractionIndex);
	m_cavityDept = cavityDept;

	m_actualBases.first = new Polygon(BASE_VERTEX_NUM);
	m_actualBases.last  = new Polygon(BASE_VERTEX_NUM);

	SetFacetParams();

	SetBaseFacets();
	CopyFacet(m_originBases.first.arr, *m_actualBases.first);
	CopyFacet(m_originBases.last.arr, *m_actualBases.last);

	SetOriginCavityPoints();

	SetCavities(m_defaultCavities);
	SetSideFacets();

	SetDefaultNormals();
	SetActualNormals();
	SetCenters();
}

void ConcaveHexagonal::RotateCavities()
{
	CavityPoints cavPoints;
	RotatePoint(m_defaultCavities.top, cavPoints.top);
	RotatePoint(m_defaultCavities.bottom, cavPoints.bottom);
	SetCavities(cavPoints);
}

void ConcaveHexagonal::Rotate(double beta, double gamma, double alpha)
{
	RotateHexagonal(beta, gamma, alpha);
	RotateCavities();
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
	SetSideFacetParams(BASE_VERTEX_NUM, 2*BASE_VERTEX_NUM);

	facetNum = BASE_VERTEX_NUM + BASE_FACET_NUM * BASE_VERTEX_NUM;

	// top facet (triangles)
	for (int i = 0; i < BASE_VERTEX_NUM; ++i)
	{
		facets[i].polygon.size = CAVITY_FACET_VERTEX_NUM;
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

void ConcaveHexagonal::SetCavities(const CavityPoints &cavities)
{
	Point3f *baseTop = m_actualBases.first->arr;
	Point3f *baseBottom = m_actualBases.last->arr;

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
	m_defaultCavities.top	 = Point3f(0, 0,  m_halfHeight - m_cavityDept);
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
