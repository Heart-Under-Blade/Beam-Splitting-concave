#include "Hexagonal.h"

const int Hexagonal::BASE_FACET_NUM;
const int Hexagonal::SIDE_VERTEX_NUMBER;
const int Hexagonal::BASE_VERTEX_NUM;

Hexagonal::Hexagonal() {}

Hexagonal::Hexagonal(double radius, double halfHeight, const complex &refractionIndex)
	: Particle(radius, halfHeight, refractionIndex)
{
	SetFacetParams();

	SetOriginBaseFacets();

	// set original positions of the base facets
	CopyFacet(m_originBases.top, facets[0]);
	CopyFacet(m_originBases.bottom, facets[7]);

	SetSideFacets(facets[0].polygon, facets[7].polygon, 1, facetNum-1);

	SetOriginNormals();
}

void Hexagonal::SetOriginNormals()
{
	double cos30 = sqrt(3)/2;

	// base facets
	m_originNormals[0] = Point3f(0, 0, -1);
	m_originNormals[7] = Point3f(0, 0, 1);

	// side facets
	m_originNormals[1] = Point3f(-cos30, -0.5, 0);
	m_originNormals[2] = Point3f(0, -1, 0);
	m_originNormals[3] = Point3f(cos30, -0.5, 0);
	m_originNormals[4] = Point3f(cos30, 0.5, 0);
	m_originNormals[5] = Point3f(0, 1, 0);
	m_originNormals[6] = Point3f(-cos30, 0.5, 0);

	// current normals
	CopyPoints(m_originNormals, normals, facetNum);

	SetDParams();
	SetExternalNormals();
}

void Hexagonal::SetFacetParams()
{
	facetNum = BASE_VERTEX_NUM + BASE_FACET_NUM;

	// base facet number
	facets[0].size = BASE_VERTEX_NUM;
	facets[facetNum-1].size = BASE_VERTEX_NUM;

	// side facet number
	for (int i = 1; i < facetNum-1; ++i)
	{
		facets[i].size = SIDE_VERTEX_NUMBER;
	}
}

void Hexagonal::SetOriginBaseFacets()
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
	facet = m_originBases.top;
	SetTwoDiagonalPoints(0, halfRadius, incircleRadius, halfHeight);
	SetTwoDiagonalPoints(1, -halfRadius, incircleRadius, halfHeight);
	SetTwoDiagonalPoints(2, -m_radius, 0, halfHeight);

	// bottom base facet
	facet = m_originBases.bottom;
	SetTwoDiagonalPoints(0, m_radius, 0, -halfHeight);
	SetTwoDiagonalPoints(1, halfRadius, -incircleRadius, -halfHeight);
	SetTwoDiagonalPoints(2, -halfRadius, -incircleRadius, -halfHeight);
}

void Hexagonal::SetSideFacets(Point3f *baseTop, Point3f *baseBottom,
							  int startIndex, int endIndex)
{
	Point3f *facet;
	int endPointIndex = BASE_VERTEX_NUM-1;

	int i1 = endPointIndex;
	int i2 = 0;

	for (int i = startIndex; i < endIndex; ++i)
	{
		facet = facets[i].polygon;

		facet[0] = baseTop[i2];
		facet[1] = baseTop[i1];

		facet[2] = baseBottom[endPointIndex-i1];
		facet[3] = baseBottom[endPointIndex-i2];

		i1 = i2;
		++i2;
	}
}

void Hexagonal::RotateBaseFacets(Point3f *baseTop, Point3f *baseBottom)
{
	for (int i = 0; i < BASE_VERTEX_NUM; ++i)
	{
		RotatePoint(m_originBases.top[i], baseTop[i]);
		RotatePoint(m_originBases.bottom[i], baseBottom[i]);
	}
}

void Hexagonal::Rotate(double beta, double gamma, double alpha)
{
	SetRotateMatrix(beta, gamma, alpha);
	RotateBaseFacets(facets[0].polygon, facets[7].polygon);
	SetSideFacets(facets[0].polygon, facets[7].polygon, 1, facetNum-1);
	RotateNormals();
}
