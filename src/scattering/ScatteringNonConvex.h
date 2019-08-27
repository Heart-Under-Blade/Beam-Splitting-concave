#pragma once

#include "Scattering.h"

/** NOTE: пучки выходят со случайно ориентированным порядком вершин */
class ScatteringNonConvex : public Scattering
{
public:
	ScatteringNonConvex(Particle *particle, Light *incidentLight,
						bool isOpticalPath, int nActs);

	void ScatterLight(double beta, double gamma, std::vector<Beam> &scaterredBeams) override;
	void ScatterLight(double beta, double gamma, const std::vector<std::vector<int>> &tracks,
							 std::vector<Beam> &scaterredBeams) override;

	double MesureOpticalPath(const Beam &beam, const Point3f sourcePoint, const std::vector<int> &track) override;
private:
	void SortFacets_faster(const Point3f &beamDir, IntArray &facetIDs);
	int FindClosestVertex(const Polygon &facet, const Point3f &beamDir);
	void CutBeamByFacet(const Facet &facet, Beam &beam,
						PolygonArray &result);

	double CalcMinDistanceToFacet(const Polygon &polygon, const Point3f &beamDir);
	void SortFacets(const Point3f &beamDir, IntArray &facetIds); ///< use 'Fast sort' algorithm

	void CutExternalBeam(const Beam &beam, std::vector<Beam> &scaterredBeams);

	void FindVisibleFacets(const Beam &beam, IntArray &facetIds);
	void FindVisibleFacetsForLight(IntArray &facetIDs);

	void SelectVisibleFacets(const Beam &beam, IntArray &facetIDs);
	void SelectVisibleFacetsForLight(IntArray &facetIDs);

	bool SetOpticalBeamParams(const Facet &facet, const Beam &incidentBeam,
							  Beam &inBeam, Beam &outBeam);

	void IntersectWithFacet(const IntArray &facetIds, int prevFacetNum,
							PolygonArray &resFacets);

	void SplitLightToBeams();

	bool IsOutgoingBeam(Beam &incidentBeam);

	int FindFacetId(int facetId, const IntArray &arr);

	void TraceFirstBeamFixedFacet(int facetID, bool &isIncident);

	void PushBeamsToTree(int facetId, const PolygonArray &polygons,
						 Beam &inBeam, Beam &outBeam);

	bool IsVisibleFacet(int facetID, const Beam &beam);

	void SplitByFacet(const IntArray &facetIDs, int facetIndex);

	bool SplitBeamByFacet(const Polygon &intersection, int facetId,
						  Beam &beam);

	void PushBeamsToBuffer(int facetID, const Beam &beam, bool hasOutBeam,
						   Beam &inBeam, Beam &outBeam, std::vector<Beam> &passed);

	void CutPolygonByFacets(const Polygon &pol,
							const IntArray &facetIds, size_t size,
							const Vector3f &polNormal, const Vector3f &clipNormal,
							const Vector3f &dir,
							PolygonArray &pols);

	void PushBeamPartsToTree(const Beam &beam,
							 const PolygonArray &parts);
	template<class T>
	void PushBeamToTree(Beam &beam, const Beam &oldBeam,
						const T &newId, int facetId,
						Location loc);

protected:
	void SplitBeams(std::vector<Beam> &scaterredBeams);

private:
	PolygonArray m_polygonBuffer;
};

