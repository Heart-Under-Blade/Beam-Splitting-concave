#pragma once

#include "Tracing.h"
#include "clipper.hpp"

struct BeamInfoConcave : public BeamInfo
{
	bool isExternal;

	BeamInfoConcave() {}
	BeamInfoConcave(Beam b, int fi, int d, bool ie)
	{
		beam = b;
		facetIndex = fi;
		dept = d;
		isExternal = ie;
	}
};

class TracingConcave : public Tracing
{
public:
	TracingConcave(Particle *particle, const Point3f &startBeamDir,
				   bool isOpticalPath, const Point3f &polarizationBasis,
				   int interReflectionNumber);

	void SplitBeamByParticle(std::vector<Beam> &outBeams,
							 double &lightSurfaceSquare) override;

	void SplitBeamByParticle(const std::vector<std::vector<int>> &tracks,
							 std::vector<Beam> &outBeams) override;

private:
	double MeasureMinDistanceToFacet(int facetIndex, const Point3f &beamDir);
	void SortFacets(int number, const Point3f &beamDir, int *facetIndices); ///< use fast sort algorithm
	void SelectVisibleFacets(const BeamInfoConcave &info, int *facetIndices,
							 int &indicesNumber);
	void SetBeamShapesByClipping(int *facetIndices, int previewFacetCount, Beam &inBeam, Beam &outBeam);
	void SetPolygonByFacet(int facetIndex, ClipperLib::Paths &polygon);
	void CutShadowsOutOfFacet(int *facetIndices, int currIndex, ClipperLib::Paths &result);

	void ProjectPointToFacet(const Point3f &point, const Point3f &direction,
							 const Point3f &facetNormal, Point3f &projection);
	void ProjectFacetToFacet(int a_index, const Point3f &normal,
							 ClipperLib::Paths &projection);

protected:
	void TraceInternalReflections(BeamInfo *tree, int size,
								  std::vector<Beam> &outBeams) override;
};
