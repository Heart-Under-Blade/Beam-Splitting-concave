#pragma once

#include "Particle.h"

class CertainAggregate : public Particle
{
public:
	CertainAggregate(const complex &refrIndex, double sizeIndex);

	void GetParticalFacetIdRangeByFacetId(int id, int &begin, int &end) const override;

protected:
	void SetFacetParams() override;
};
