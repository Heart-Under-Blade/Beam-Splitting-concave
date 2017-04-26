#include "Hexagonal.h"
#include "global.h"

Hexagonal::Hexagonal() {}

Hexagonal::Hexagonal(double radius, double halfHeight, const complex &refractionIndex)
	: Particle(radius, halfHeight, refractionIndex)
{
	m_symmetryAngle = M_PI/3;

	SetFacetParams();
	SetBaseFacets(defaultState.facets[0], defaultState.facets[7]);

	SetSideFacets(defaultState.facets[0], defaultState.facets[7]);
	SetDefaultNormals();
	SetActualState();
}

void Hexagonal::SetFacetParams()
{
	facetNum = BASE_VERTEX_NUM + BASE_FACET_NUM;
	SetSideFacetParams(1, facetNum-1);

	// base facet number
	defaultState.facets[0].polygon.size = BASE_VERTEX_NUM;
	defaultState.facets[facetNum-1].polygon.size = BASE_VERTEX_NUM;
}

void Hexagonal::SetSideFacetParams(int first, int last)
{
	m_sideFacetIDs = Couple<int>{first, last};

	for (int i = first; i < last; ++i)
	{
		defaultState.facets[i].polygon.size = SIDE_VERTEX_NUM;
	}
}

void Hexagonal::SetBaseFacets(Facet &baseTop, Facet &baseBottom)
{
	Point3f *facet;

	int halfNumber = BASE_VERTEX_NUM/2;
	double halfRadius = m_radius/2;
	double incircleRadius = (sqrt(3) * m_radius) / 2;

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
	facet = baseTop.polygon.arr;
	SetTwoDiagonalPoints(0, halfRadius, incircleRadius, m_halfHeight);
	SetTwoDiagonalPoints(1,-halfRadius, incircleRadius, m_halfHeight);
	SetTwoDiagonalPoints(2,-m_radius, 0,  m_halfHeight);

	// bottom base facet
	facet = baseBottom.polygon.arr;
	SetTwoDiagonalPoints(0, m_radius, 0, -m_halfHeight);
	SetTwoDiagonalPoints(1, halfRadius, -incircleRadius, -m_halfHeight);
	SetTwoDiagonalPoints(2,-halfRadius, -incircleRadius, -m_halfHeight);
}

void Hexagonal::SetSideFacets(Facet &baseTop, Facet &baseBottom)
{
	const Point3f *top = baseTop.polygon.arr;
	const Point3f *bot = baseBottom.polygon.arr;

	int endPointIndex = BASE_VERTEX_NUM-1;

	int i1 = endPointIndex;
	int i2 = 0;

	for (int i = m_sideFacetIDs.first; i <  m_sideFacetIDs.last; ++i)
	{
		Point3f *facet = defaultState.facets[i].polygon.arr;

		facet[0] = top[i2];
		facet[1] = top[i1];

		facet[2] = bot[endPointIndex-i1];
		facet[3] = bot[endPointIndex-i2];

		i1 = i2;
		++i2;
	}
}
