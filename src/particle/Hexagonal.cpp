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
	CopyPoints(m_originBases.top, facets[0], vertexNums[0]);
	CopyPoints(m_originBases.bottom, facets[7], vertexNums[7]);

	SetSideFacets(facets[0], facets[7], 1, facetNum-1);

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
	vertexNums[0] = BASE_VERTEX_NUM;
	vertexNums[facetNum-1] = BASE_VERTEX_NUM;

	// side facet number
	for (int i = 1; i < facetNum-1; ++i)
	{
		vertexNums[i] = SIDE_VERTEX_NUMBER;
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
		facet = facets[i];

		facet[0] = baseTop[i2];
		facet[1] = baseTop[i1];

		facet[2] = baseBottom[endPointIndex-i1];
		facet[3] = baseBottom[endPointIndex-i2];

		i1 = i2;
		++i2;
	}
}

/// OPT: добавить потом
//void Hexagon::Rotate(double beta, double gamma)
//{
//#ifdef __FASTMATH_H
//	double cosB, cosG,
//			sinB, sinG;

//	sincos(beta, &sinB, &cosB);
//	sincos(gamma, &sinG, &cosG);
//#else
//	double cosB = cos(beta);
//	double cosG = cos(gamma);

//	double sinB = sin(beta);
//	double sinG = sin(gamma);
//#endif

//	/// OPT: дооптимизировать
//	double m11 = cosB*cosG;
//	double m21 = sinG;
//	double m31 = -sinB*cosG;

//	double m12 = -cosB*sinG;
//	double m22 = cosG;
//	double m32 = sinB*sinG;

//	double m13 = sinB;
//	double m23 = 0;
//	double m33 = cosB;

//	for (int i = 0; i < m_facetNum; ++i)
//	{
//		for (int j = 0; j < m_vertexNums[i]; ++j)
//		{
//			Point3 oldVertex = m_facets[i][j];
//			m_facets[i][j].X = oldVertex.X*m11 + oldVertex.Y*m12 + oldVertex.Z*m13;
//			m_facets[i][j].Y = oldVertex.X*m21 + oldVertex.Y*m22 + oldVertex.Z*m23;
//			m_facets[i][j].Z = oldVertex.X*m31 + oldVertex.Y*m32 + oldVertex.Z*m33;
//		}
//	}
//}

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
	RotateBaseFacets(facets[0], facets[7]);
	SetSideFacets(facets[0], facets[7], 1, facetNum-1);
	RotateNormals();
}
