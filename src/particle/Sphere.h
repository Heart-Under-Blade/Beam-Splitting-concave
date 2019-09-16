#pragma once

#include "Particle.h"

class Sphere : public Particle
{
public:
	Sphere(const complex &refrIndex, int nParallels, int nMeridians,
		   double radius);
};
