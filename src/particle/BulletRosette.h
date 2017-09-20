#pragma once

//#include "Particle.h"
#include "Hexagonal.h"

class BulletRosette : public Hexagonal
{
public:
	BulletRosette();
	BulletRosette(const complex &refrIndex, double diameter, double height,
				  double peakHeight);

	// Particle interface
public:
	void GetAggPartFacetIDRange(int id, int &begin, int &end) const override;
};
