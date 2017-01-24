#pragma once

#include "Tracing.h"
#include "BeamClipper.h"

/** NOTE: пучки выходят со случайно ориентированным порядком вершин */
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
	void CutIncidentBeam2(int facetId, Beam &beam, bool &isDivided);
	double CalcMinDistanceToFacet(const Polygon &polygon, const Point3f &beamDir);
	void SortFacets(const Point3f &beamDir, IntArray &facetIds); ///< use 'Fast sort' algorithm

	void CutShadowsFromFacet2(int facetId, const IntArray &facetIds, int prevFacetNum,
							 const Beam &beam, Polygon *resFacets, int &resSize);

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

	void CatchExternalBeam(const Beam &beam, std::vector<Beam> &scatteredBeams);

	void PushBeamToTree(Beam &beam, int facetId, int level, Location location);

	void FindVisibleFacets_initial(IntArray &facetIds);
	void FindVisibleFacets(const Beam &beam, IntArray &facetIds);

	void SelectVisibleFacets(const Beam &beam, IntArray &facetIds);

	void SetOpticalBeamParams(int facetId, Beam &incidentBeam,
							  Beam &inBeam, Beam &outBeam, bool &hasOutBeam);

	void IntersectWithFacet(const IntArray &facetIds, int prevFacetNum,
							Polygon *resFacets, int &resSize, bool &hasIntersection);

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

