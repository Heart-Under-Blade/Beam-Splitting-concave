#pragma once

//#include "Particle.h"
#include "Hexagonal.h"

class BulletRosette : public Hexagonal
{
public:
	BulletRosette();
	BulletRosette(const complex &refrIndex, double diameter, double height,
				  double peakHeight);
};
