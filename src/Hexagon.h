#pragma once

#include "Particle.h"

#define SIDE_NUMBER			6	///< number of side facets
#define SIDE_VERTEX_NUMBER	4	///< number of vertex of the each side facets
#define BASE_FACET_NUMBER	2	///< number of bases of hexagon

/**
 * @brief The Hexagon class
 * The prism particle with 6 number of side facets.
 */
class Hexagon : public Particle
{
public:
	Hexagon(double radius, double halfHeight, const complex &refractionIndex);

	void Rotate(double beta, double gamma, double alpha) override;

private:
	void SetBases();
	void SetSides();
	void SetFacetParams();
	void SetNormals();
	void SetExternalNormals();
	void SetDParams();
};

