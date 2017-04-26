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
	Hexagonal(double m_radius, double m_halfHeight, const complex &m_refractionIndex);

//	void Rotate(double beta, double gamma, double alpha) override;

protected:
	static const int BASE_FACET_NUM  = 2;	///< number of bases of hexagon
	static const int SIDE_VERTEX_NUM = 4;	///< number of vertex of the each side facet
	static const int BASE_VERTEX_NUM = 6;	///< number of side facets

	Couple<int> m_sideFacetIDs;

//	Couple<Polygon> m_originBases;
//	Couple<Polygon*> m_actualBases; ///< actual coordinates of base facets

protected:
//	void SetDefaultNormals() override;
	void SetFacetParams() override;

	void SetBaseFacets(Facet &baseTop, Facet &baseBottom);
	void SetSideFacets(Facet &baseTop, Facet &baseBottom);
//	void RotateHexagonal(double gamma, double beta, double alpha);
//	void SetSideNormals(int beginId);
	void SetSideFacetParams(int first, int last);

private:
//	void RotateBaseFacets();
};

