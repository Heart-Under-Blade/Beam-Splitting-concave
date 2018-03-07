#pragma once

#include "Scattering.h"

/** NOTE: пучки выходят со случайно ориентированным порядком вершин */
class ScatteringNonConvex : public Scattering
{
public:
	ScatteringNonConvex(Particle *particle, Light *incidentLight,
						bool isOpticalPath, int numberOfActs);

	void ScatterLight(double beta, double gamma, std::vector<Beam> &scaterredBeams) override;
	void ScatterLight(double beta, double gamma, const std::vector<std::vector<int>> &tracks,
							 std::vector<Beam> &scaterredBeams) override;
private:
	void SortFacets_faster(const Point3f &beamDir, IntArray &facetIDs);
	int FindClosestVertex(const Polygon &facet, const Point3f &beamDir);
	void CutBeamByFacet(int facetID, Beam &beam, bool &isDivided,
						Polygon *resultBeams, int &resultSize);

	double CalcMinDistanceToFacet(const Polygon &polygon, const Point3f &beamDir);
	void SortFacets(const Point3f &beamDir, IntArray &facetIds); ///< use 'Fast sort' algorithm

	void CutFacetByShadows(int facetID, const IntArray &shadowFacetIDs, int prevFacetNum,
						   PolygonArray &resFacets);

	void ProjectPointToFacet(const Point3f &point, const Point3f &direction,
							 const Point3f &facetNormal, Point3f &projection);

	void CatchExternalBeam(const Beam &beam, std::vector<Beam> &scatteredBeams);

	void FindVisibleFacets(const Beam &beam, IntArray &facetIds);
	void FindVisibleFacetsForLight(IntArray &facetIDs);

	void SelectVisibleFacets(const Beam &beam, IntArray &facetIDs);
	void SelectVisibleFacetsForLight(IntArray &facetIDs);

	void SetOpticalBeamParams(int facetID, const Beam &incidentBeam,
							  Beam &inBeam, Beam &outBeam, bool &hasOutBeam);

	void IntersectWithFacet(const IntArray &facetIDs, int prevFacetNum,
							PolygonArray &resFacets);

	void SplitLightToBeams();

	bool isExternalNonEmptyBeam(Beam &incidentBeam);

	int FindFacetID(int facetID, const IntArray &arr);

	void TraceFirstBeamFixedFacet(int facetID, bool &isIncident);

#ifdef _TRACK_ALLOW
//	void AddToTrack(Beam &beam, int facetId);
#endif

	void PushBeamsToTree(const Beam &beam, int facetID, bool hasOutBeam,
						 Beam &inBeam, Beam &outBeam);
	void PushBeamsToTree(int facetID, const PolygonArray &polygons,
						 Beam &inBeam, Beam &outBeam);

	bool IsVisibleFacet(int facetID, const Beam &beam);

	void SplitByFacet(const IntArray &facetIDs, int facetIndex);

	void TraceSecondaryBeamByFacet(Beam &beam, int facetID, bool &isDivided);

	void PushBeamsToBuffer(int facetID, const Beam &beam, bool hasOutBeam, Beam &inBeam, Beam &outBeam, std::vector<Beam> &passed);

protected:
	void TraceBeams(std::vector<Beam> &scaterredBeams);
};

