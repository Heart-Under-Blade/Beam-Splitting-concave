#include "Hexagonal.h"
#include "global.h"
#include <algorithm>

Hexagonal::Hexagonal() {}

Hexagonal::Hexagonal(const complex &refrIndex, double diameter, double height)
{
	SetSize(diameter, height);

	double size = std::max(m_height, m_diameter);
	Init(8, refrIndex, M_PI/3, size);

	SetFacetParams();
	SetBases(defaultState.facets[0], defaultState.facets[7]);
	SetSides(defaultState.facets[0], defaultState.facets[7]);
	SetDefaultNormals();
	SetActualState();
}

void Hexagonal::SetSize(double diameter, double height)
{
	m_diameter = diameter;
	m_height = height;
}

void Hexagonal::SetFacetParams()
{
	SetSideFacetParams(1, m_facetNum-1);

	// base facet number
	defaultState.facets[0].size = BASE_VERTEX_NUM;
	defaultState.facets[m_facetNum-1].size = BASE_VERTEX_NUM;
}

void Hexagonal::SetSideFacetParams(int first, int last)
{
	m_sideFacetIDs = Couple<int>{first, last};

	for (int i = first; i < last; ++i)
	{
		defaultState.facets[i].size = SIDE_VERTEX_NUM;
	}
}

void Hexagonal::SetBases(Facet &baseTop, Facet &baseBottom)
{
	Point3f *facet;

	double radius = m_diameter/2;
	double halfHeight = m_height/2;

	int halfNumber = BASE_VERTEX_NUM/2;
	double halfRadius = radius/2;
	double inRadius = (sqrt(3) * radius) / 2;

	auto SetTwoDiagonalPoints = [&] (int startPointIndex, double x, double y, double z)
	{
		int endPointIndex = startPointIndex + halfNumber;

		for (int i = startPointIndex; i <= endPointIndex; i += halfNumber)
		{
			facet[i].cx = x;
			facet[i].cy = y;
			facet[i].cz = z;

			x = -x;
			y = -y;
		}
	};

	// top base facet
	facet = baseTop.arr;
	SetTwoDiagonalPoints(0, halfRadius, inRadius, halfHeight);
	SetTwoDiagonalPoints(1, -halfRadius, inRadius, halfHeight);
	SetTwoDiagonalPoints(2, -radius, 0, halfHeight);

	// bottom base facet
	facet = baseBottom.arr;
	SetTwoDiagonalPoints(0, radius, 0, -halfHeight);
	SetTwoDiagonalPoints(1, halfRadius, -inRadius, -halfHeight);
	SetTwoDiagonalPoints(2, -halfRadius, -inRadius, -halfHeight);
}

void Hexagonal::SetSides(Facet &baseTop, Facet &baseBottom)
{
	const Point3f *top = baseTop.arr;
	const Point3f *bot = baseBottom.arr;

	int endPointIndex = BASE_VERTEX_NUM-1;

	int i1 = endPointIndex;
	int i2 = 0;

	for (int i = m_sideFacetIDs.first; i <  m_sideFacetIDs.last; ++i)
	{
		Point3f *facet = defaultState.facets[i].arr;

		facet[0] = top[i2];
		facet[1] = top[i1];

		facet[2] = bot[endPointIndex-i1];
		facet[3] = bot[endPointIndex-i2];

		i1 = i2;
		++i2;
	}
}
