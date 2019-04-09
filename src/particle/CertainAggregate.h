#pragma once

#include "Particle.h"

class CertainAggregate : public Particle
{
public:
	CertainAggregate();
	void GetPartByFacet(Facet *facet, Array<Facet*> &facets) override;

protected:
	void SetFacetParams() override;
};
