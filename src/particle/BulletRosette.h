#pragma once

//#include "Particle.h"
#include "Column.h"

class BulletRosette : public Particle
{
public:
	BulletRosette();
	BulletRosette(const complex &refrIndex, const Size &size,
				  double peakHeight);

	void GetParticalFacetIdRange(Facet *facet, int &begin, int &end) const override;
};
