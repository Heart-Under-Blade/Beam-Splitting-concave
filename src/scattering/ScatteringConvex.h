#pragma once

#include "Scattering.h"

class ScatteringConvex : public Scattering
{
public:
	ScatteringConvex(Particle *particle, const Light &incidentLight, int maxActNo);
        
protected:
	void SplitOriginalBeam(std::vector<Beam> &externalBeams) override;
	void SelectVisibleFacets(const Beam &beam, Array<Facet *> &facets) override;
	void PushBeamsToBuffer(Facet *facet, BeamPair<Beam> &beams,
						   std::vector<Beam> &scatteredBeams);
};
