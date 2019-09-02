#pragma once

//#include "Particle.h"
#include "Column.h"

class BulletRosette : public Particle
{
public:
	BulletRosette();
	BulletRosette(const Size &size, double peakHeight);

	void GetPartByFacet(Facet *facet, Array<Facet *> &facets) override;
};
