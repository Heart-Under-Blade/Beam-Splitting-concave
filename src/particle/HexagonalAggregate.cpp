#include "HexagonalAggregate.h"
#include "common.h"

HexagonalAggregate::HexagonalAggregate(const complex &refrIndex, const Size &size,
									   int particleNumber)
	: Column(8*particleNumber, refrIndex, size, true)
{
	SetSymmetry(M_PI, 2*M_PI);
	SetFacetParams();

//  SetBases1 REF: вынести в отд. ф-цию
	{
		Facet &baseTop = elems[0].origin;
		Facet &baseBottom = elems[7].origin;

		Point3f *facet;

		double radius = m_size.diameter/2;
		double halfHeight = m_size.height/2;

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
		Facet &baseTop = elems[8].origin;
		Facet &baseBottom = elems[15].origin;

		Point3f *facet;

		double radius = m_size.diameter/2;
		double halfHeight = m_size.height/2;

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
	SetSides(elems[0].origin, elems[7].origin);

	SetSideFacetParams(9, 15);
	SetSides(elems[8].origin, elems[15].origin);

	SetDefaultNormals();
	SetDefaultCenters();
	Reset();
}

void HexagonalAggregate::SetFacetParams()
{
	elems[0].origin.nVertices = BASE_VERTEX_NUM;
	elems[7].origin.nVertices = BASE_VERTEX_NUM;

	elems[8].origin.nVertices = BASE_VERTEX_NUM;
	elems[15].origin.nVertices = BASE_VERTEX_NUM;

	elems[1].actual.isOverlayedOut = true;
	elems[2].actual.isOverlayedOut = true;
	elems[3].actual.isOverlayedOut = true;
	elems[7].actual.isOverlayedOut = true;
	elems[10].actual.isOverlayedOut = true;
	elems[11].actual.isOverlayedOut = true;
	elems[12].actual.isOverlayedOut = true;
	elems[15].actual.isOverlayedOut = true;

	for (int i = 0; i < nElems; ++i) // OPT: кол-во затеняемых гарней на самом деле меньше
	{
		elems[i].actual.isOverlayedIn = true;
	}
}
