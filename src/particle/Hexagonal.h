#pragma once

#include "Particle.h"

struct ParticleParams {
	double radius;
	double halfHeight;
	complex refractionIndex;
};

/**
 * @brief The Hexagon class
 * The prism particle with 6 number of side facets.
 */
class Hexagonal : public Particle
{
public:
	Hexagonal(const ParticleParams &params);
	Hexagonal(double m_radius, double halfHeight, const complex &refractionIndex);

	void Rotate(double beta, double gamma, double alpha) override;

protected:
	static const int BASE_FACET_NUMBER = 2;		///< number of bases of hexagon
	static const int SIDE_VERTEX_NUMBER = 4;	///< number of vertex of the each side facet
	static const int BASE_VERTEX_NUMBER = 6;	///< number of side facets

	struct Bases
	{
		Point3f top[BASE_VERTEX_NUMBER];
		Point3f bottom[BASE_VERTEX_NUMBER];
	}
	m_originBases;

protected:
	void SetOriginNormals() override;
	void SetFacetParams() override;

	void SetOriginBaseFacets();
	void SetSideFacets(Point3f *baseTop, Point3f *baseBottom,
					   int startIndex, int endIndex);
	void RotateBaseFacets(Point3f *baseTop, Point3f *baseBottom);
};

