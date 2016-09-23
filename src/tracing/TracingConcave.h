#pragma once

#include "Tracing.h"
#include "clipper.hpp"

class TracingConcave : public Tracing
{
public:
	TracingConcave(Particle *particle, const Point3f &startBeamDir,
				   bool isOpticalPath, const Point3f &polarizationBasis,
				   int interReflectionNumber);

	void SplitBeamByParticle(std::vector<OutBeam> &outBeams, double &lightSurfaceSquare) override;

private:
	double m_startBeamDParam;

private:
	double MeasureMinDistanceToFacet(int facetIndex);
	void SortFacets(int *facetIndices, int number); ///< use fast sort algorithm
	void SelectVisibleFacets(int *facetIndices, int &indicesNumber);
	void SetBeamShapesByClipping(int *facetIndices, int previewFacetCount, Beam &inBeam, Beam &outBeam);
	void SetPolygonByFacet(int facetIndex, ClipperLib::Paths &polygon);
	void CutShadesOutOfFacet(int *facetIndices, int previewFacetCount, ClipperLib::Paths &result);
};
