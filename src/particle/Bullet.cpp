#include "Bullet.h"
#include "global.h"

Bullet::Bullet()
{

}

Bullet::Bullet(const complex &refrIndex, double diameter, double height, double peakHeight)
{
	SetSize(diameter, height);

	double size = std::max(m_height + peakHeight, m_diameter);
	Init(13, refrIndex, size);

	SetSymmetry(M_PI/2, M_PI/3);
	SetFacetParams();

	Facet baseTop;
	SetBases(baseTop, defaultFacets[12]);
	SetSides(baseTop, defaultFacets[12]);
	Point3f peak = Point3f(0, 0, m_height/2 + peakHeight);
	SetPeakFacets(0, 6, baseTop.arr, peak); // top facets (triangles)

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
						   const Point3f &cavityPoint)
{
	// base facet point indices
	int p0 = BASE_VERTEX_NUM-1;
	int p1 = 0;

	for (int i = start; i < end; ++i)
	{
		defaultFacets[i].arr[0] = baseFacet[p0];
		defaultFacets[i].arr[1] = baseFacet[p1];
		defaultFacets[i].arr[2] = cavityPoint;
		p0 = p1;
		++p1;
	}
}

void Bullet::SetFacetParams()
{
	for (int i = 0; i < facetNum; ++i)
	{
		facets[i].isVisibleIn = false;
		facets[i].isVisibleOut = false;
	}

	for (int i = 0; i < BASE_VERTEX_NUM; ++i)
	{
		defaultFacets[i].size = 3;
	}

	SetSideFacetParams(6, 12);

	defaultFacets[12].size = 6;
}
