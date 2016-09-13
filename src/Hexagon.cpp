#include "Hexagon.h"

#include "vector_lib.h"

Hexagon::Hexagon(double radius, double halfHeight, const complex &refractionIndex)
	: Particle(radius, halfHeight, refractionIndex)
{
	SetFacetParams();
	SetBases();
	SetSides();
	SetNormals();
}

void Hexagon::SetNormals()
{
	double cos30 = sqrt(3)/2;

	originNormals[0] = normals[0] = Point3f(0, 0, -1);
	originNormals[1] = normals[1] = Point3f(-cos30, -0.5, 0);
	originNormals[2] = normals[2] = Point3f(0, -1, 0);
	originNormals[3] = normals[3] = Point3f(cos30, -0.5, 0);
	originNormals[4] = normals[4] = Point3f(cos30, 0.5, 0);
	originNormals[5] = normals[5] = Point3f(0, 1, 0);
	originNormals[6] = normals[6] = Point3f(-cos30, 0.5, 0);
	originNormals[7] = normals[7] = Point3f(0, 0, 1);

	SetDParams();
	SetExternalNormals();
}

void Hexagon::SetDParams()
{
	for (int i = 0; i < facetNum; ++i)
	{
		double d = DotProduct(facets[i][0], normals[i]);
		normals[i].D_PARAM = -d;
	}
}

void Hexagon::SetExternalNormals()
{
	for (int i = 0; i < facetNum; ++i)
	{
		externalNormals[i].X = -normals[i].X;
		externalNormals[i].Y = -normals[i].Y;
		externalNormals[i].Z = -normals[i].Z;
		externalNormals[i].D_PARAM = -normals[i].D_PARAM;
	}
}

void Hexagon::SetFacetParams()
{
	facetNum = SIDE_NUMBER + BASE_FACET_NUMBER;

	// bases
	vertexNums[0] = SIDE_NUMBER;
	vertexNums[facetNum-1] = SIDE_NUMBER;

	// sides
	for (int i = 1; i < facetNum-1; ++i)
	{
		vertexNums[i] = SIDE_VERTEX_NUMBER;
	}
}

void Hexagon::SetBases()
{
	Point3f *facet;
	Point3f *originFacet;

	int halfRadius = radius/2;
	int halfPointNumber = SIDE_NUMBER/2;
	int incircleRadius = (sqrt(3) * radius) / 2;

	auto SetTwoDiagonalPoints = [&] (int startPointIndex, int x, int y, int z)
	{
		int endPointIndex = startPointIndex + halfPointNumber;

		for (int i = startPointIndex; i <= endPointIndex; i += halfPointNumber)
		{
			originFacet[i].X = facet[i].X = x;
			originFacet[i].Y = facet[i].Y = y;
			originFacet[i].Z = facet[i].Z = z;

			x *= -1;
			y *= -1;
		}
	};

	// 0 - base facet
	facet = facets[0];
	originFacet = originFacets[0];
	SetTwoDiagonalPoints(0, halfRadius, incircleRadius, halfHeight);
	SetTwoDiagonalPoints(1, -halfRadius, incircleRadius, halfHeight);
	SetTwoDiagonalPoints(2, -radius, 0, halfHeight);

	// 7 - base facet
	facet = facets[7];
	originFacet = originFacets[7];
	SetTwoDiagonalPoints(0, radius, 0, -halfHeight);
	SetTwoDiagonalPoints(1, halfRadius, -incircleRadius, -halfHeight);
	SetTwoDiagonalPoints(2, -halfRadius, -incircleRadius, -halfHeight);
}

void Hexagon::SetSides()
{
	Point3f *facet;
	Point3f *originFacet;

	Point3f *baseTop = facets[0];
	Point3f *baseBottom = facets[7];

	int endPointIndex = SIDE_NUMBER-1;

	int i1 = endPointIndex;
	int i2 = 0;

	for (int i = 1; i <= SIDE_NUMBER; ++i)
	{
		facet = facets[i];
		originFacet = originFacets[i];

		originFacet[0] = facet[0] = baseTop[i2];
		originFacet[1] = facet[1] = baseTop[i1];

		originFacet[2] = facet[2] = baseBottom[endPointIndex-i1];
		originFacet[3] = facet[3] = baseBottom[endPointIndex-i2];

		i1 = i2;
		++i2;
	}
}

/// TODO: добавить потом
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

//	/// TODO: дооптимизировать
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

void Hexagon::Rotate(double beta, double gamma, double alpha)
{
#ifdef __FASTMATH_H
	double cosA, cosB, cosG,
			sinA, sinB, sinG;

	sincos(alpha, &sinA, &cosA);
	sincos(beta, &sinB, &cosB);
	sincos(gamma, &sinG, &cosG);
#else
	double cosA = cos(alpha);
	double cosB = cos(beta);
	double cosG = cos(gamma);

	double sinA = sin(alpha);
	double sinB = sin(beta);
	double sinG = sin(gamma);
#endif

	/// TODO: дооптимизировать
	double cosAcosB = cosA*cosB;

	double m11 = cosAcosB*cosG - sinA*sinG;
	double m21 = sinA*cosB*cosG + cosA*sinG;
	double m31 = -sinB*cosG;

	double m12 = -(cosAcosB*sinG + sinA*cosG);
	double m22 = cosA*cosG - sinA*cosB*sinG;
	double m32 = sinB*sinG;

	double m13 = cosA*sinB;
	double m23 = sinA*sinB;
	double m33 = cosB;

	auto RotatePoint = [&] (const Point3f &point, Point3f &result)
	{
		result.X = point.X*m11 + point.Y*m12 + point.Z*m13;
		result.Y = point.X*m21 + point.Y*m22 + point.Z*m23;
		result.Z = point.X*m31 + point.Y*m32 + point.Z*m33;
	};

	auto RotateFacet = [&] (int facetIndex)
	{
		for (int j = 0; j < vertexNums[facetIndex]; ++j)
		{
			RotatePoint(originFacets[facetIndex][j], facets[facetIndex][j]);
		}
	};

	// resolve vertices of base facets
	RotateFacet(0);
	RotateFacet(7);

	// resolve normals
	for (int i = 0; i < facetNum; ++i)
	{
		RotatePoint(originNormals[i], normals[i]);
	}

	SetSides();
	SetDParams();
	SetExternalNormals();
}
