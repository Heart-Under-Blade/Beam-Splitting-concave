#include "Hexagonal.h"
#include "global.h"

Hexagonal::Hexagonal() {}

Hexagonal::Hexagonal(double radius, double halfHeight, const complex &refractionIndex)
	: Particle(radius, halfHeight, refractionIndex)
{
	m_symmetryAngle = M_PI/3;
//	m_actualBases.first = &(facets[0].polygon);
//	m_actualBases.last  = &(facets[7].polygon);

	SetFacetParams();
	SetBaseFacets(defaultState.facets[0], defaultState.facets[7]);

	// set original positions of the base facets
//	CopyFacet(m_originBases.first.arr, facets[0]);
//	CopyFacet(m_originBases.last.arr, facets[7]);

	SetSideFacets(defaultState.facets[0], defaultState.facets[7]);
	SetDefaultNormals();
	SetActualState();
}

//void Hexagonal::SetSideNormals(int beginId)
//{
//	double cos30 = sqrt(3)/2;
//	m_originNormals[beginId++] = Point3f(-cos30,-0.5, 0);
//	m_originNormals[beginId++] = Point3f(0,-1, 0);
//	m_originNormals[beginId++] = Point3f(cos30,-0.5, 0);
//	m_originNormals[beginId++] = Point3f(cos30, 0.5, 0);
//	m_originNormals[beginId++] = Point3f(0, 1, 0);
//	m_originNormals[beginId++] = Point3f(-cos30, 0.5, 0);
//}

//void Hexagonal::SetDefaultNormals()
//{
//	// base facets
//	m_originNormals[0] = Point3f(0, 0,-1);
//	m_originNormals[7] = Point3f(0, 0, 1);

//	SetSideNormals(1);
//}

void Hexagonal::SetFacetParams()
{
	SetSideFacetParams(1, facetNum-1);
	facetNum = BASE_VERTEX_NUM + BASE_FACET_NUM;

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

//void Hexagonal::RotateBaseFacets()
//{
//	for (int i = 0; i < BASE_VERTEX_NUM; ++i)
//	{
//		RotatePoint(m_originBases.first.arr[i], m_actualBases.first->arr[i]);
//		RotatePoint(m_originBases.last.arr[i], m_actualBases.last->arr[i]);
//	}
//}

//void Hexagonal::RotateHexagonal(double beta, double gamma, double alpha)
//{
//	Particle::Rotate(beta, gamma, alpha);
//	RotateBaseFacets();
//	SetSideFacets();
//}

//void Hexagonal::Rotate(double beta, double gamma, double alpha)
//{
//	RotateHexagonal(beta, gamma, alpha);
//	RotateNormals();
//}
