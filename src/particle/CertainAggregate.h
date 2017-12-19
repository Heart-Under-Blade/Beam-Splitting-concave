#pragma once

#include "Particle.h"

class CertainAggregate : public Particle
{
public:
	CertainAggregate(const complex &refrIndex, double sizeIndex);

	void GetAggPartFacetIDRange(int id, int &begin, int &end) const override;

protected:
	void SetFacetParams() override;

	// Particle interface
private:
	void Resize(double sizeIndex);
};
