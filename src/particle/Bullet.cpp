#include "Bullet.h"
#include "common.h"

Bullet::Bullet()
{
}

Bullet::Bullet(const Size &size, double peakHeight)
	: Column(13, size, false)
{
	SetSymmetry(M_PI/2, M_PI/3);
	SetFacetParams();

	Facet baseTop;
	SetBases(baseTop, elems[7].original);
	SetSides(baseTop, elems[7].original);
	Point3f peak = Point3f(0, 0, m_size.height/2 + peakHeight);
	SetPeakFacets(8, 13, baseTop.vertices, peak); // top facets (triangles)

	SetDefaultNormals();
	SetDefaultCenters();
	ResetPosition();
}

void Bullet::SetBaseFacet(Facet &facet)
{
	double radius = m_size.diameter/2;
	double halfHeight = m_size.height/2;

	double halfRadius = radius/2;
	double inRadius = (sqrt(3) * radius) / 2;

	// top base facet
	Point3f *polygon = facet.vertices;
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
	elems[0].original.vertices[0] = baseFacet[p0];
	elems[0].original.vertices[1] = baseFacet[p1];
	elems[0].original.vertices[2] = peakPoint;
	p0 = p1;
	++p1;

	for (int i = start; i != end; ++i)
	{
		elems[i].original.vertices[0] = baseFacet[p0];
		elems[i].original.vertices[1] = baseFacet[p1];
		elems[i].original.vertices[2] = peakPoint;
		p0 = p1;
		++p1;
	}
}

void Bullet::SetFacetParams()
{
	for (int i = 0; i < nElems; ++i)
	{
		elems[i].actual.isOverlayedIn = true;
		elems[i].actual.isOverlayedOut = true;
	}

	elems[0].original.nVertices = 3;

	for (int i = 8; i < 13; ++i)
	{
		elems[i].original.nVertices = 3;
	}

	SetSideFacetParams(1, 7);

	elems[7].original.nVertices = BASE_VERTEX_NUM;
}
