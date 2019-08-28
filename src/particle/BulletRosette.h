#pragma once

//#include "Particle.h"
#include "Hexagonal.h"

class BulletRosette : public Particle
{
public:
	BulletRosette();
	BulletRosette(const complex &refrIndex, double diameter, double height,
				  double peakHeight);

	void GetParticalFacetIdRangeByFacetId(int id, int &begin, int &end) const override;
};
