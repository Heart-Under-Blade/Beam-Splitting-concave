#pragma once

#include "Tracing.h"
#include "clipper.hpp"

typedef std::vector<std::vector<Point3f>> Polygons;

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

	double BeamCrossSection(const Beam &beam) const override;

	void SplitBeamByParticle(std::vector<Beam> &outBeams,
							 double &lightSurfaceSquare) override;

	void SplitBeamByParticle(const std::vector<std::vector<int>> &tracks,
							 std::vector<Beam> &outBeams) override;

private:
	double MeasureMinDistanceToFacet(int facetIndex, const Point3f &beamDir);
	void SortFacets(int number, const Point3f &beamDir, int *facetIndices); ///< use fast sort algorithm
	void SelectVisibleFacets(const BeamInfoConcave &info, int *facetIndices,
							 int &indicesNumber);
	void SetBeamShapesByClipping(int *facetIndices, int facetCount, bool isExternal, Beam &inBeam, Beam &outBeam);
	void SetPolygonByFacet(const Point3f *facet, int size, ClipperLib::Paths &polygon);
	void CutShadowsOutOfFacet(int *facetIndices, int facetCount, const Point3f &normal,
							  ClipperLib::Paths &result);

	void ProjectPointToFacet(const Point3f &point, const Point3f &direction,
							 const Point3f &facetNormal, Point3f &projection);
	void ProjectFacetToFacet(int a_index, const Point3f &normal,
							 ClipperLib::Paths &projection);

	void SetBeamShapeByPolygon(Beam &beam, const ClipperLib::Paths &result);

	void CutBeamShape(const Beam &outBeam, Beam &incidentBeam);

	void DivideConcavePolygon(const Point3f *polygon, int size,
							  const Point3f &normal) const;

	void findDividePoint(const Point3f *polygon, int size,
						 int i0, int i1, const Point3f &normal,
						 Point3f &x, int &nextPointIndex) const;

	void FillSubpolygon(int begin, int end,
						const Point3f *polygon, int size,
						std::vector<Point3f> &subpolygon) const;

protected:
	void TraceInternalReflections(BeamInfoConcave *tree, int size,
								  std::vector<Beam> &outBeams);
};
