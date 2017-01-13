#pragma once

#include "Tracing.h"
#include "BeamClipper.h"

typedef std::vector<std::vector<Point3f>> Polygons;

/** NOTE: пучки выходят со случайно ориентированными вершинами */
class TracingConcave : public Tracing
{
public:
	TracingConcave(Particle *particle, const Point3f &startBeamDir,
				   bool isOpticalPath, const Point3f &polarizationBasis,
				   int interReflectionNumber);

	double BeamCrossSection(const Beam &beam) const override;

	void SplitBeamByParticle(std::vector<Beam> &outBeams) override;

	void SplitBeamByParticle(const std::vector<std::vector<int>> &tracks,
							 std::vector<Beam> &outBeams) override;
private:
	BeamClipper m_clipper;

private:
	void CutIncidentBeam(int facetId, Beam &beam, bool &isDivided);
	void CutIncidentBeam2(int facetId, Beam &beam, bool &isDivided);
	double CalcMinDistanceToFacet(const Polygon &polygon, const Point3f &beamDir);
	void SortFacets(const Point3f &beamDir, IntArray &facetIds); ///< use fast sort algorithm
	void CutShadowsFromFacet(int facetId, IntArray facetIds, int handledFacetNum,
							 const Beam &beam, ClipperLib::Paths &resultFacet);
	void CutShadowsFromFacet2(int facetId, IntArray facetIds, int handledFacetNum,
							 const Beam &beam, Polygon *allFacets, int &allSize);

	void ProjectPointToFacet(const Point3d &point, const Point3d &direction,
							 const Point3d &facetNormal, Point3d &projection);

	void ProjectPointToFacet(const Point3f &point, const Point3f &direction,
							 const Point3f &facetNormal, Point3f &projection);

	void ProjectFacetToFacet(const Polygon &a_facet, const Point3f &a_dir,
							 const Point3f &b_normal,
							 ClipperLib::Path &projection);


	void CutBeamByFacet(ClipperLib::Paths &beamPolygon, int facetId,
						const Point3f &direction, const Point3f &facetNormal,
						ClipperLib::Paths &result);

	// recursive
	void DividePolygon(const std::vector<Point3f> &polygon,
					   const Point3f &normal, Polygons &polygons) const;
	void DivideConcavePolygon(const Point3f *polygon, int size,
							  const Point3f &normal,
							  Polygons &polygons) const;
	void FindDividePoint(const std::vector<Point3f> &polygon,
						 int i0, int i1, const Point3f &normal,
						 Point3f &x, int &nextPointIndex) const;
	void FillSubpolygon(int begin, int end,
						const std::vector<Point3f> &polygon,
						std::vector<Point3f> &subpolygon) const;

//	bool isOrderReversed(const Point3f oldNormal, const ClipperLib::Path polygon);
//	void InversePolygonOrder(ClipperLib::Path &polygon);

	void CatchExternalBeam(const Beam &beam, std::vector<Beam> &scatteredBeams);

	void PushBeamToTree(Beam &beam, int facetId, int level, Location location);

	void FindVisibleFacets_initial(IntArray &facetIds);
	void FindVisibleFacets(const Beam &beam, IntArray &facetIds);

	void SelectVisibleFacets(const Beam &beam, IntArray &facetIds);

	void SetOpticalBeamParams(int facetId, Beam &incidentBeam,
							  Beam &inBeam, Beam &outBeam, bool &hasOutBeam);

	void IntersectWithFacet(const IntArray &facetIds, int handledFacetNum,
								Polygon *allFacets, int &allSize, bool &hasIntersection);

	void TraceFirstBeam();

	void SelectFirstBeamVisibleFacets(IntArray &facetIds);

	bool HasExternalBeam(Beam &incidentBeam);

#ifdef _TRACK_ALLOW
	void PrintTrack(const Beam &beam, int facetId);
	void AddToTrack(Beam &beam, int facetId);
#endif

	void PushBeamsToTree(int level, int facetId, bool hasOutBeam, Beam &inBeam, Beam &outBeam);

protected:
	void TraceSecondaryBeams(std::vector<Beam> &scaterredBeams);
};

