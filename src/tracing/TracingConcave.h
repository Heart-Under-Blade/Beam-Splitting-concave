#pragma once

#include "Tracing.h"

/** NOTE: пучки выходят со случайно ориентированным порядком вершин */
class TracingConcave : public Tracing
{
public:
	TracingConcave(Particle *particle, const Point3f &startBeamDir,
				   bool isOpticalPath, const Point3f &polarizationBasis,
				   int interReflectionNumber);

	double BeamCrossSection(const Beam &beam) const override;

	void SplitBeamByParticle(std::vector<Beam> &scaterredBeams) override;
	void SplitBeamByParticle(const std::vector<std::vector<int>> &tracks,
							 std::vector<Beam> &scaterredBeams) override;

private:
	void CutBeamByFacet(int facetId, Beam &beam, bool &isDivided);
	double CalcMinDistanceToFacet(const Polygon &polygon, const Point3f &beamDir);
	void SortFacets(const Point3f &beamDir, IntArray &facetIds); ///< use 'Fast sort' algorithm

	void CutFacetByShadows(int facetID, const IntArray &shadowFacetIDs, int prevFacetNum,
						   PolygonArray &resFacets);

	void ProjectPointToFacet(const Point3f &point, const Point3f &direction,
							 const Point3f &facetNormal, Point3f &projection);

	void CatchExternalBeam(const Beam &beam, std::vector<Beam> &scatteredBeams);

	void PushBeamToTree(Beam &beam, int facetId, int level, Location location);
	void PushBeamToTree(Beam &beam);

	void FindVisibleFacetsForWavefront(IntArray &facetIds);
	void FindVisibleFacets2(const Beam &beam, IntArray &facetIds,
							bool isWavefront = false);
	void FindVisibleFacets(const Beam &beam, IntArray &facetIds);

	void SelectVisibleFacets(const Beam &beam, IntArray &facetIds);

	void SetOpticalBeamParams(int facetId, Beam &incidentBeam,
							  Beam &inBeam, Beam &outBeam, bool &hasOutBeam);

	void IntersectWithFacet(const IntArray &facetIDs, int prevFacetNum,
							PolygonArray &resFacets);

	void TraceFirstBeam();

	void SelectVisibleFacetsForWavefront(IntArray &facetIds);

	bool HasExternalBeam(Beam &incidentBeam);

#ifdef _TRACK_ALLOW
	void PrintTrack(const Beam &beam, int facetId);
	void AddToTrack(Beam &beam, int facetId);
#endif

	void PushBeamsToTree(int level, int facetID, bool hasOutBeam,
						 Beam &inBeam, Beam &outBeam);

	bool IsVisibleFacet(int facetID, const Beam &beam);

	void PushBeamsToTree(int facetID, const PolygonArray &polygons,
						 Beam &inBeam, Beam &outBeam);

protected:
	void TraceSecondaryBeams(std::vector<Beam> &scaterredBeams);
};

