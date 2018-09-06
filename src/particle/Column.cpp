#include "Column.h"
#include "global.h"
#include <algorithm>

Column::Column() {}

Column::Column(size_t nFacets, const complex &refrIndex, const Size &size,
			   bool isNonConvex)
	: Particle(nFacets, refrIndex, isNonConvex),
	  m_size(size)
{
}

void Column::SetFacetParams()
{
	SetSideFacetParams(1, nFacets-1);

	// base facet number
	defaultFacets[0].nVertices = BASE_VERTEX_NUM;
	defaultFacets[nFacets-1].nVertices = BASE_VERTEX_NUM;
}

void Column::SetSideFacetParams(int first, int last)
{
	m_sideFacetIDs = Couple<int>{first, last};// REF: как-нибудь избавиться от этого

	for (int i = first; i < last; ++i)
	{
		defaultFacets[i].nVertices = SIDE_VERTEX_NUM;
	}
}

void Column::SetBases(Facet &top, Facet &bottom)
{
	Point3f *facet;

	double radius = m_size.diameter/2;
	double halfHeight = m_size.height/2;

	double halfRadius = radius/2;
	double inRadius = (sqrt(3) * radius) / 2;

	// top base facet
	facet = top.arr;
	SetTwoDiagonalPoints(0, facet, halfRadius, inRadius, halfHeight);
	SetTwoDiagonalPoints(1, facet, -halfRadius, inRadius, halfHeight);
	SetTwoDiagonalPoints(2, facet, -radius, 0, halfHeight);

	// bottom base facet
	facet = bottom.arr;
	SetTwoDiagonalPoints(0, facet, radius, 0, -halfHeight);
	SetTwoDiagonalPoints(1, facet, halfRadius, -inRadius, -halfHeight);
	SetTwoDiagonalPoints(2, facet, -halfRadius, -inRadius, -halfHeight);
}

void Column::SetTwoDiagonalPoints(int index, Point3f *facet,
									 double x, double y, double z)
{
	int halfNumber = BASE_VERTEX_NUM/2;
	int endIndex = index + halfNumber;
	facet[index] = Point3f(x, y, z);
	facet[endIndex] = Point3f(-x, -y, z);
}

void Column::SetSides(Facet &baseTop, Facet &baseBottom)
{
	const Point3f *top = baseTop.arr;
	const Point3f *bot = baseBottom.arr;

	int endIndex = BASE_VERTEX_NUM-1;

	int i1 = endIndex;
	int i2 = 0;

	for (int i = m_sideFacetIDs.first; i < m_sideFacetIDs.last; ++i)
	{
		Point3f *facet = defaultFacets[i].arr;

		facet[0] = top[i2];
		facet[1] = top[i1];

		facet[2] = bot[endIndex-i1];
		facet[3] = bot[endIndex-i2];

		i1 = i2;
		++i2;
	}
}
