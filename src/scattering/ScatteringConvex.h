#pragma once

#include "Scattering.h"

class ScatteringConvex : public Scattering
{
public:
	ScatteringConvex(Particle *particle, const Light &incidentLight,
					 int maxActNo, const complex &m_refractiveIndex);

	void ScatterLight(TrackNode *trackTree, std::vector<Beam> &scatteredBeams);

	double MeasureOpticalPath(const Beam &beam, const Point3f sourcePoint,
							 const std::vector<int> &track) override;

	double MeasureFullOpticalPath(const Beam &beam, const Point3f sourcePoint,
								  const std::vector<int> &track) override;
protected:
	void SplitOriginalBeam(std::vector<Beam> &externalBeams) override;
	void SelectVisibleFacets(const Beam &beam, Array<Facet *> &facets) override;
	void PushBeamsToBuffer(Facet *facet, BeamPair<Beam> &beams,
						   std::vector<Beam> &scatteredBeams);
private:
	void SetWorkFacetsFromTracks();
};
