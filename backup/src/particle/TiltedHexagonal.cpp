//#include "TiltedHexagonal.h"
//#include "global.h"
//#include <iostream>

//const int TiltedHexagonal::BASE_FACET_NUM;
//const int TiltedHexagonal::SIDE_VERTEX_NUMBER;
//const int TiltedHexagonal::BASE_VERTEX_NUM;

//TiltedHexagonal::TiltedHexagonal(double radius, double halfHeight, const complex &refractionIndex, double tiltAngle)
////	: Particle(radius, halfHeight, refractionIndex)
//{
//	m_symmetryAngle = 2*M_PI;
//	m_tiltAngle = tiltAngle;

//	SetFacetParams();
//	SetBaseFacets();
//	TiltBaseFacets();

//	// set original positions of the base facets
//	CopyFacet(m_originBases[0].arr, facets[0]);
//	CopyFacet(m_originBases[1].arr, facets[7]);

//	SetSideFacets(facets[0].polygon.arr, facets[7].polygon.arr, 1, facetNum-1);
////	SetDefaultNormals();
//	SetActualNormals();
//}

//void TiltedHexagonal::TiltBaseFacets()
//{
//	double tilt = 15.0*M_PI/180.0,
//		   tan_ang = tan(m_tiltAngle*M_PI/180.0);

//	double h[6];
//	h[0] = m_radius * cos(5.0*M_PI/3.0-tilt) * tan_ang;
//	h[1] = m_radius * cos(tilt) * tan_ang,
//	h[2] = m_radius * cos(M_PI/3.0-tilt) * tan_ang,
//	h[3] = m_radius * cos(2.0*M_PI/3.0-tilt) * tan_ang,
//	h[4] = m_radius * cos(M_PI-tilt) * tan_ang,
//	h[5] = m_radius * cos(4.0*M_PI/3.0-tilt) * tan_ang;

//	int endPointIndex = BASE_VERTEX_NUM-1;

//	for (int i = 0; i < facetNum; ++i)
//	{
//		m_originBases[0].arr[i] += h[i];
//		m_originBases[1].arr[endPointIndex-i] += h[i];
//	}
//}

////void TiltedHexagonal::SetSideNormals(int beginId)
////{
////	double cos30 = sqrt(3)/2;
////	m_originNormals[beginId++] = Point3f(-cos30,-0.5, 0);
////	m_originNormals[beginId++] = Point3f(0,-1, 0);
////	m_originNormals[beginId++] = Point3f(cos30,-0.5, 0);
////	m_originNormals[beginId++] = Point3f(cos30, 0.5, 0);
////	m_originNormals[beginId++] = Point3f(0, 1, 0);
////	m_originNormals[beginId++] = Point3f(-cos30, 0.5, 0);
////}

////void TiltedHexagonal::SetDefaultNormals()
////{
//	// base facets
////	m_originNormals[0] = -m_originBases[0].Normal();
////	m_originNormals[7] = -m_originBases[1].Normal();

////	SetSideNormals(1);
////}

//void TiltedHexagonal::SetFacetParams()
//{
//	facetNum = BASE_VERTEX_NUM + BASE_FACET_NUM;

//	// base facet number
//	facets[0].polygon.size = BASE_VERTEX_NUM;
//	facets[facetNum-1].polygon.size = BASE_VERTEX_NUM;

//	// side facet number
//	for (int i = 1; i < facetNum-1; ++i)
//	{
//		facets[i].polygon.size = SIDE_VERTEX_NUMBER;
//	}
//}

//void TiltedHexagonal::SetBaseFacets()
//{
//	Point3f *facet;

//	int halfNumber = BASE_VERTEX_NUM/2;
//	double halfRadius = m_radius/2;
//	double incircleRadius = (sqrt(3) * m_radius) / 2;

//	auto SetTwoDiagonalPoints = [&] (int startPointIndex, double x, double y, double z)
//	{
//		int endPointIndex = startPointIndex + halfNumber;

//		for (int i = startPointIndex; i <= endPointIndex; i += halfNumber)
//		{
//			facet[i].cx = x;
//			facet[i].cy = y;
//			facet[i].cz = z;

//			x = -x;
//			y = -y;
//		}
//	};

//	// top base facet
//	facet = m_originBases[0].arr;
//	SetTwoDiagonalPoints(0, halfRadius, incircleRadius, m_halfHeight);
//	SetTwoDiagonalPoints(1,-halfRadius, incircleRadius, m_halfHeight);
//	SetTwoDiagonalPoints(2,-m_radius, 0,  m_halfHeight);

//	// bottom base facet
//	facet = m_originBases[1].arr;
//	SetTwoDiagonalPoints(0, m_radius, 0, -m_halfHeight);
//	SetTwoDiagonalPoints(1, halfRadius, -incircleRadius, -m_halfHeight);
////	std::cout << "!!111" << facet[2].cx << std::endl // DEB
////			  << facet[5].cx << std::endl;
//	SetTwoDiagonalPoints(2,-halfRadius, -incircleRadius, -m_halfHeight);
//}

//void TiltedHexagonal::SetSideFacets(Point3f *baseTop, Point3f *baseBottom,
//							  int startIndex, int endIndex)
//{
//	int endPointIndex = BASE_VERTEX_NUM-1;

//	int i1 = endPointIndex;
//	int i2 = 0;

//	for (int i = startIndex; i < endIndex; ++i)
//	{
//		Point3f *facet = facets[i].polygon.arr;

//		facet[0] = baseTop[i2];
//		facet[1] = baseTop[i1];

//		facet[2] = baseBottom[endPointIndex-i1];
//		facet[3] = baseBottom[endPointIndex-i2];

//		i1 = i2;
//		++i2;
//	}
//}

////void TiltedHexagonal::RotateBaseFacets(Point3f *baseTop, Point3f *baseBottom)
////{
////	for (int i = 0; i < BASE_VERTEX_NUM; ++i)
////	{
////		RotatePoint(m_originBases[0].arr[i], baseTop[i]);
////		RotatePoint(m_originBases[1].arr[i], baseBottom[i]);
////	}
////}

////void TiltedHexagonal::Rotate(double beta, double gamma, double alpha)
////{
////	Particle::Rotate(beta, gamma, alpha);

////	Point3f *baseTop	= facets[0].polygon.arr;
////	Point3f *baseBottom = facets[7].polygon.arr;

////	RotateBaseFacets(baseTop, baseBottom);
////	SetSideFacets(baseTop, baseBottom, 1, facetNum-1);

////	RotateNormals();
////}
