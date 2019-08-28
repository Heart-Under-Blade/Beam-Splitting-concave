#include "Bullet.h"
#include "global.h"

Bullet::Bullet()
{

}

Bullet::Bullet(const complex &refrIndex, double diameter, double height, double peakHeight)
{
	isConcave = false;
	SetSize(diameter, height);
	Init(13, refrIndex);

	SetSymmetry(M_PI/2, M_PI/3);
	SetFacetParams();

	Facet baseTop;
	SetBases(baseTop, defaultFacets[7]);
	SetSides(baseTop, defaultFacets[7]);
	Point3f peak = Point3f(0, 0, m_height/2 + peakHeight);
	SetPeakFacets(8, 13, baseTop.arr, peak); // top facets (triangles)

	SetDefaultNormals();
	SetDefaultCenters();
	Reset();
}

void Bullet::SetBaseFacet(Facet &facet)
{
	double radius = m_diameter/2;
	double halfHeight = m_height/2;

	double halfRadius = radius/2;
	double inRadius = (sqrt(3) * radius) / 2;

	// top base facet
	Point3f *polygon = facet.arr;
	SetTwoDiagonalPoints(0, polygon, halfRadius, inRadius, halfHeight);
	SetTwoDiagonalPoints(1, polygon, -halfRadius, inRadius, halfHeight);
	SetTwoDiagonalPoints(2, polygon, -radius, 0, halfHeight);
}

void Bullet::SetPeakFacets(int start, int end, const Point3f *baseFacet,
						   const Point3f &peakPoint)
{
	// base facet point indices
	int p0 = BASE_VERTEX_NUM-1;
	int p1 = 0;
	defaultFacets[0].arr[0] = baseFacet[p0];
	defaultFacets[0].arr[1] = baseFacet[p1];
	defaultFacets[0].arr[2] = peakPoint;
	p0 = p1;
	++p1;

	for (int i = start; i != end; ++i)
	{
		defaultFacets[i].arr[0] = baseFacet[p0];
		defaultFacets[i].arr[1] = baseFacet[p1];
		defaultFacets[i].arr[2] = peakPoint;
		p0 = p1;
		++p1;
	}
}

void Bullet::SetFacetParams()
{
	for (int i = 0; i < nFacets; ++i)
	{
		facets[i].isVisibleIn = false;
		facets[i].isVisibleOut = false;
	}

	defaultFacets[0].nVertices = 3;

	for (int i = 8; i < 13; ++i)
	{
		defaultFacets[i].nVertices = 3;
	}

	SetSideFacetParams(1, 7);

	defaultFacets[7].nVertices = 6;
}
