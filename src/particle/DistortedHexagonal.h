#pragma once

#include "Hexagonal.h"

/**
 * @brief The Hexagon class
 * The prism particle with 6 number of side facets.
 */
class DistortedHexagonal : public Hexagonal
{
public:
	DistortedHexagonal();
	DistortedHexagonal(const complex &refrIndex, double diameter, double height,
					   double angle);

private:
	void DistortBases(double angle);
};

