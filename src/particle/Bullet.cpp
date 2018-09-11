#include "Bullet.h"
#include "global.h"

Bullet::Bullet()
{

}

Bullet::Bullet(const complex &refrIndex, const Size &size, double peakHeight)
	: Column(13, refrIndex, size, false)
{
	SetSymmetry(M_PI/2, M_PI/3);
	SetFacetParams();

	Facet baseTop;
	SetBases(baseTop, elems[7].origin);
	SetSides(baseTop, elems[7].origin);
	Point3f peak = Point3f(0, 0, m_size.height/2 + peakHeight);
	SetPeakFacets(8, 13, baseTop.arr, peak); // top facets (triangles)

	SetDefaultNormals();
	SetDefaultCenters();
	Reset();
}

void Bullet::SetBaseFacet(Facet &facet)
{
	double radius = m_size.diameter/2;
	double halfHeight = m_size.height/2;

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
	elems[0].origin.arr[0] = baseFacet[p0];
	elems[0].origin.arr[1] = baseFacet[p1];
	elems[0].origin.arr[2] = peakPoint;
	p0 = p1;
	++p1;

	for (int i = start; i != end; ++i)
	{
		elems[i].origin.arr[0] = baseFacet[p0];
		elems[i].origin.arr[1] = baseFacet[p1];
		elems[i].origin.arr[2] = peakPoint;
		p0 = p1;
		++p1;
	}
}

void Bullet::SetFacetParams()
{
	for (int i = 0; i < nElems; ++i)
	{
		elems[i].actual.isOverlayedIn = false;
		elems[i].actual.isOverlayedOut = false;
	}

	elems[0].origin.nVertices = 3;

	for (int i = 8; i < 13; ++i)
	{
		elems[i].origin.nVertices = 3;
	}

	SetSideFacetParams(1, 7);

	elems[7].origin.nVertices = BASE_VERTEX_NUM;
}
