#pragma once

#include "Scattering.h"

/** NOTE: пучки выходят со случайно ориентированным порядком вершин */
class ScatteringNonConvex : public Scattering
{
public:
	ScatteringNonConvex(Particle *particle, const Light &incidentLight,
						int maxActNo, const complex &refractiveIndex);
	~ScatteringNonConvex();

	void ScatterLight(TrackNode *trackTree, std::vector<Beam> &scatteredBeams);

protected:
	void SplitOriginalBeam(std::vector<Beam> &scatteredBeams) override;

	void ReleaseBeam(Beam &beam) override;
	bool IsFinalAct(const Beam &beam) override;
	bool isFinalFacet(int index, Array<Facet*> &facets) override;
	void PushBeamsToBuffer(Beam &parentBeam, Facet *facet,
						   bool hasOutBeam) override;
	void SelectVisibleFacets(const Beam &beam, Array<Facet*> &visibleFacets) override;

	double MeasureOpticalPath(const Beam &beam, const Point3f sourcePoint,
							  const std::vector<int> &track) override;

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
	/**
	 * @brief Sorts facets of the Particle by distance from the beginnig of
	 * latest beam splitting in descending order.
	 * @param beamDir beam direction
	 * @param facets unordered facets
	 */
	void SortFacetsByDistance(const Vector3f &beamDir, Array<Facet*> &facets) const;

	/**
	 * @brief Finds vertex of the Facet polygon that located more near than
	 * others to beginnig of the Beam.
	 * @param facet facet polygon
	 * @param beamDir beam direction
	 * @return index of the closest vertex in the beam polygon
	 */
	int FindClosestVertex(const Polygon &facet, const Point3f &beamDir) const;

	bool FindRestOfBeamShape(Facet *facet, const Beam &beam, PolygonStack &rest);

	double CalcMinDistanceToFacet(Polygon *facet, const Point3f &beamDir);
	void SortFacets(const Point3f &beamDir, Array<Facet*> &facets); ///< use 'Fast sort' algorithm

	bool FindVisiblePartOfFacet(const Array<Facet*> &facets, int nCheckedFacets,
								 PolygonStack &pols);

	void PushBeamsToTree(Facet *facet, BeamPair<Beam> &beams,
						 const PolygonStack &polygons,
						 std::vector<Beam> &scatteredBeams);

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
						bool loc);

protected:
	void SplitBeams(std::vector<Beam> &scaterredBeams);

private:
	PolygonArray m_polygonBuffer;
};

