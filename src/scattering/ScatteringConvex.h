#pragma once

#include "Scattering.h"

class ScatteringConvex : public Scattering
{
public:
	ScatteringConvex(const complex &refractiveIndex,
					 int maxActNo, double minBeamEnergy = MIN_BEAM_NRG,
					 double farFresnelZone = FAR_FRESNEL_ZONE);

	void Scatter(TrackNode *trackTree, std::vector<Beam> &scatteredBeams);

	double MeasureOpticalPath(const BeamInfo &info,
							  const Point3f sourcePoint) override;

	double MeasureFullOpticalPath(const BeamInfo &info,
								  const Point3f sourcePoint) override;
protected:
	void SplitStartBeam() override;
	void SplitSecondaryBeams() override;
	void SelectVisibleFacets(const Beam &beam, Array<Facet *> *facets) override;
	void ResolveExternalBeam(const Beam &beam) override;

private:
	void SetWorkFacetsFromTracks();

	Point3f lastPoint;
};
