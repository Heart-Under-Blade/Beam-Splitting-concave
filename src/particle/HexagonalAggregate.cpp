#include "HexagonalAggregate.h"
#include "global.h"

HexagonalAggregate::HexagonalAggregate(const complex &refrIndex, double diameter, double height,
									   int particleNumber)
{
	isConcave = true;
	SetSize(diameter, height);
	Init(8*particleNumber, refrIndex);

	SetSymmetry(M_PI, 2*M_PI);
	SetFacetParams();

//  SetBases1 REF: вынести в отд. ф-цию
	{
		Facet &baseTop = defaultFacets[0];
		Facet &baseBottom = defaultFacets[7];

		Point3f *facet;

		double radius = m_diameter/2;
		double halfHeight = m_height/2;

		double offset = halfHeight + radius + 1/*REF: вынести как конст.*/;

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
		facet = baseTop.arr;
		SetTwoDiagonalPoints(0, halfRadius, inRadius, halfHeight + offset);
		SetTwoDiagonalPoints(1, -halfRadius, inRadius, halfHeight + offset);
		SetTwoDiagonalPoints(2, -radius, 0, halfHeight + offset);

		// bottom base facet
		facet = baseBottom.arr;
		SetTwoDiagonalPoints(0, radius, 0, -halfHeight + offset);
		SetTwoDiagonalPoints(1, halfRadius, -inRadius, -halfHeight + offset);
		SetTwoDiagonalPoints(2, -halfRadius, -inRadius, -halfHeight + offset);
	}

//	SetBases2 REF: вынести в отд. ф-цию
	{
		Facet &baseTop = defaultFacets[8];
		Facet &baseBottom = defaultFacets[15];

		Point3f *facet;

		double radius = m_diameter/2;
		double halfHeight = m_height/2;

		double offset = halfHeight + radius + 1/*PART_OFFSET*/;

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
				z = -z;
			}
		};

		// top base facet
		facet = baseTop.arr;
		SetTwoDiagonalPoints(2, halfRadius, halfHeight+offset, inRadius);
		SetTwoDiagonalPoints(1, -halfRadius, halfHeight+offset, inRadius);
		SetTwoDiagonalPoints(0, -radius, halfHeight+offset, 0);

		// bottom base facet
		facet = baseBottom.arr;
		SetTwoDiagonalPoints(2, radius, -halfHeight+offset, 0);
		SetTwoDiagonalPoints(1, halfRadius, -halfHeight+offset, -inRadius);
		SetTwoDiagonalPoints(0, -halfRadius, -halfHeight+offset, -inRadius);
	}

	SetSideFacetParams(1, 7);
	SetSides(defaultFacets[0], defaultFacets[7]);

	SetSideFacetParams(9, 15);
	SetSides(defaultFacets[8], defaultFacets[15]);

	SetDefaultNormals();
	SetDefaultCenters();
	Reset();
}

void HexagonalAggregate::SetFacetParams()
{
	defaultFacets[0].nVertices = BASE_VERTEX_NUM;
	defaultFacets[7].nVertices = BASE_VERTEX_NUM;

	defaultFacets[8].nVertices = BASE_VERTEX_NUM;
	defaultFacets[15].nVertices = BASE_VERTEX_NUM;

	facets[1].isVisibleOut = false;
	facets[2].isVisibleOut = false;
	facets[3].isVisibleOut = false;
	facets[7].isVisibleOut = false;
	facets[10].isVisibleOut = false;
	facets[11].isVisibleOut = false;
	facets[12].isVisibleOut = false;
	facets[15].isVisibleOut = false;

	for (int i = 0; i < nFacets; ++i) // OPT: кол-во затеняемых гарней на самом деле меньше
	{
		facets[i].isVisibleIn = false;
	}
}
