#pragma once

#include "Scattering.h"

class ScatteringConvex : public Scattering
{
public:
	ScatteringConvex(Particle *particle, const Light &incidentLight, int maxActNo);

protected:
	void SplitOriginBeam(std::vector<Beam> &scatteredBeams) override;
	void SelectVisibleFacets(const Beam &beam, Array<Facet *> &facets) override;
	void PushBeamsToBuffer(Facet *facet, Splitting &splitting,
						   std::vector<Beam> &scatteredBeams);
};
