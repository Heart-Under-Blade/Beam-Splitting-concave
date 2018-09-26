#pragma once

#include "Particle.h"

struct Size
{
	Size() {}
	Size(double d, double h) {diameter = d; height = h;}
	double diameter;
	double height;
};

/**
 * @brief The Column class
 */
class Column : public Particle
{
public:
	Column();
	Column(int nFacets, const complex &refrIndex, const Size &size,
		   bool isNonConvex);

protected:
	static const int BASE_NUM = 2;			///< number of bases
	static const int SIDE_VERTEX_NUM = 4;	///< number of vertex of the each side facet
	static const int BASE_VERTEX_NUM = 6;	///< number of side facets

protected:
	Couple<int> m_sideFacetIDs;
	Size m_size;

protected:
	void SetFacetParams() override;
	void SetBases(Facet &top, Facet &bottom);
	void SetTwoDiagonalPoints(int index, Point3f *facet,
							  double x, double y, double z);
	void SetSides(Facet &baseTop, Facet &baseBottom);
	void SetSideFacetParams(int first, int last);
};

