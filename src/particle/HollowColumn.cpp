#include "HollowColumn.h"
#include "common.h"
#include <algorithm>

HollowColumn::HollowColumn(const Size &size, double cavityAngle)
	: Column(18, size, true)
{
	double angleD = Orientation::DegToRad(cavityAngle);
	double r = m_size.diameter/2;
	m_cavityDept = (r*sin(angleD))/cos(angleD);

	SetSymmetry(M_PI/2, M_PI/3);
	SetFacetParams();

	Facet baseTop, baseBottom;
	SetBases(baseTop, baseBottom);

	SetOriginCavityPoints();

	SetCavities(baseTop, baseBottom, m_defaultStateCavities);
	SetSides(baseTop, baseBottom);

	SetDefaultNormals();
	ResetPosition();
	SetDefaultCenters();
}

void HollowColumn::SetFacetParams()
{
	SetSideFacetParams(BASE_VERTEX_NUM, 2*BASE_VERTEX_NUM);

	// top facet (triangles)
	for (int i = 0; i < BASE_VERTEX_NUM; ++i)
	{
		elems[i].original.nVertices = CAVITY_FACET_VERTEX_NUM;
		elems[i].actual.isOverlayedOut = true;
	}

	// bottom facet (triangles)
	for (int i = 2*BASE_VERTEX_NUM; i < 3*BASE_VERTEX_NUM; ++i)
	{
		elems[i].original.nVertices = CAVITY_FACET_VERTEX_NUM;
		elems[i].actual.isOverlayedOut = true;
	}

	for (int i = BASE_VERTEX_NUM; i < 2*BASE_VERTEX_NUM; ++i)
	{
		elems[i].actual.isOverlayedIn = true;
	}
}

void HollowColumn::SetCavities(Facet &baseTop, Facet &baseBottom,
								   const CavityPoints &cavities)
{
	Point3f *bTop = baseTop.vertices;
	Point3f *bBot = baseBottom.vertices;

	SetCavityFacets(0, BASE_VERTEX_NUM, bTop, cavities.top); // top facets (triangles)
	SetCavityFacets(2*BASE_VERTEX_NUM, 3*BASE_VERTEX_NUM, bBot, cavities.bottom); // bottom facets (triangles)
}

void HollowColumn::SetCavityFacets(int start, int end,
								   const Point3f *baseFacet,
								   const Point3f &cavityPoint)
{
	// base facet point indices
	int p0 = BASE_VERTEX_NUM-1;
	int p1 = 0;

	for (int i = start; i < end; ++i)
	{
		elems[i].original.vertices[0] = baseFacet[p0];
		elems[i].original.vertices[1] = baseFacet[p1];
		elems[i].original.vertices[2] = cavityPoint;
		p0 = p1;
		++p1;
	}
}

void HollowColumn::SetOriginCavityPoints()
{
	m_defaultStateCavities.top	  = Point3f(0, 0,  m_size.height/2 - m_cavityDept);
	m_defaultStateCavities.bottom = Point3f(0, 0, -m_size.height/2 + m_cavityDept);
}
