#pragma once

#include "Particle.h"

class CertainAggregate : public Particle
{
public:
	CertainAggregate(const complex &refrIndex, double sizeIndex);

	// Particle interface
protected:
	void SetFacetParams() override;

	// Particle interface
public:
	void GetAggPartFacetIDRange(int id, int &begin, int &end) const override;
private:
        void Resize(double sizeIndex);
};
