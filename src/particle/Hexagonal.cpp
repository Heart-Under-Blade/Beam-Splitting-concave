#include "Hexagonal.h"
#include "global.h"
#include <algorithm>

Hexagonal::Hexagonal() {}

Hexagonal::Hexagonal(const complex &refrIndex, double diameter, double height)
{
	SetSize(diameter, height);

	double size = std::max(m_height, m_diameter);
	Init(8, refrIndex, M_PI/3, M_PI/2, size);

	SetFacetParams();

	SetBases(defaultFacets[0], defaultFacets[7]);
	SetSides(defaultFacets[0], defaultFacets[7]);

	SetDefaultNormals();
	SetActualState();
	SetDefaultCenters();
}

void Hexagonal::SetSize(double diameter, double height)
{
	m_diameter = diameter;
	m_height = height;
}

void Hexagonal::SetFacetParams()
{
	SetSideFacetParams(1, facetNum-1);

	// base facet number
	defaultFacets[0].size = BASE_VERTEX_NUM;
	defaultFacets[facetNum-1].size = BASE_VERTEX_NUM;
}

void Hexagonal::SetSideFacetParams(int first, int last)
{
	m_sideFacetIDs = Couple<int>{first, last};// REF: как-нибудь избавиться от этого

	for (int i = first; i < last; ++i)
	{
		defaultFacets[i].size = SIDE_VERTEX_NUM;
	}
}

void Hexagonal::SetBases(Facet &top, Facet &bottom)
{
	Point3f *facet;

	double radius = m_diameter/2;
	double halfHeight = m_height/2;

	int halfNumber = BASE_VERTEX_NUM/2;
	double halfRadius = radius/2;
	double inRadius = (sqrt(3) * radius) / 2;

	auto SetTwoDiagonalPoints = [&] (int startIndex, double x, double y, double z)
	{
		int endIndex = startIndex + halfNumber;

		for (int i = startIndex; i <= endIndex; i += halfNumber)
		{
			facet[i] = Point3f(x, y, z);
			x = -x;
			y = -y;
		}
	};

	// top base facet
	facet = top.arr;
	SetTwoDiagonalPoints(0, halfRadius, inRadius, halfHeight);
	SetTwoDiagonalPoints(1, -halfRadius, inRadius, halfHeight);
	SetTwoDiagonalPoints(2, -radius, 0, halfHeight);

	// bottom base facet
	facet = bottom.arr;
	SetTwoDiagonalPoints(0, radius, 0, -halfHeight);
	SetTwoDiagonalPoints(1, halfRadius, -inRadius, -halfHeight);
	SetTwoDiagonalPoints(2, -halfRadius, -inRadius, -halfHeight);
}

void Hexagonal::SetSides(Facet &baseTop, Facet &baseBottom)
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
