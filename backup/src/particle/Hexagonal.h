#pragma once

#include "Particle.h"

/**
 * @brief The Hexagon class
 * The prism particle with 6 side facets.
 */
class Hexagonal : public Particle
{
public:
	Hexagonal();
	Hexagonal(const complex &refrIndex, double diameter, double height);

protected:
	static const int BASE_NUM		 = 2;	///< number hexagon's of bases
	static const int SIDE_VERTEX_NUM = 4;	///< number of vertex of the each side facet
	static const int BASE_VERTEX_NUM = 6;	///< number of side facets

protected:
	Couple<int> m_sideFacetIDs;
	double m_diameter;
	double m_height;

protected:
	void SetFacetParams() override;
	void SetBases(Facet &top, Facet &bottom);
	void SetTwoDiagonalPoints(int index, Point3f *facet,
							  double x, double y, double z);
	void SetSides(Facet &baseTop, Facet &baseBottom);
	void SetSideFacetParams(int first, int last);
	void SetSize(double diameter, double height);
};

