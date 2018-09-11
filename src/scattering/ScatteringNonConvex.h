#pragma once

#include "Scattering.h"

/** NOTE: пучки выходят со случайно ориентированным порядком вершин */
class ScatteringNonConvex : public Scattering
{
public:
	ScatteringNonConvex(Particle *particle, Light *incidentLight,
						bool isOpticalPath, int nActs);

	void ScatterLight(std::vector<Beam> &scaterredBeams) override;
	void ScatterLight(const std::vector<std::vector<int>> &tracks,
					  std::vector<Beam> &scaterredBeams) override;
private:
	void SortFacets_faster(const Point3f &beamDir, Array<Facet*> &facets);
	int FindClosestVertex(const Polygon &facet, const Point3f &beamDir);
	void CutBeamByFacet(Facet *facet, Beam &beam,
						PolygonArray &result);

	double CalcMinDistanceToFacet(Polygon *facet, const Point3f &beamDir);
	void SortFacets(const Point3f &beamDir, Array<Facet*> &facets); ///< use 'Fast sort' algorithm

	void CutExternalBeam(const Beam &beam, std::vector<Beam> &scaterredBeams);

	void FindVisibleFacets(const Beam &beam, Array<Facet *> &facets);
	void FindVisibleFacetsForLight(Array<Facet*> &facets);

	void SelectVisibleFacets(const Beam &beam, Array<Facet *> &facets);
	void SelectVisibleFacetsForLight(Array<Facet*> &facetIDs);

	bool SetOpticalBeamParams(Facet *facet, const Beam &incidentBeam,
							  Beam &inBeam, Beam &outBeam);

	void IntersectWithFacet(const Array<Facet*> &facets, size_t nCheckedFacets,
							PolygonArray &resFacets);

	void SplitLightToBeams();

	bool IsOutgoingBeam(Beam &incidentBeam);

	int FindFacetId(int facetId, const Array<Facet*> &arr);

	void TraceFirstBeamFixedFacet(int facetID, bool &isIncident);

	void PushBeamsToTree(Facet *facet, const PolygonArray &polygons,
						 Beam &inBeam, Beam &outBeam);

	bool IsVisibleFacet(Facet *facet, const Beam &beam);

	void SplitByFacet(const Array<Facet*> &facets, size_t nCheckedFacets);

	bool SplitBeamByFacet(const Polygon &intersection,
						  Facet *facet, Beam &beam);

	void PushBeamsToBuffer(int facetID, const Beam &beam, bool hasOutBeam,
						   Beam &inBeam, Beam &outBeam, std::vector<Beam> &passed);

	void CutPolygonByFacets(const Polygon &pol,
							const Array<Facet*> &facets, size_t size,
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
	void ScatterBeams(std::vector<Beam> &scaterredBeams);

private:
	PolygonArray m_polygonBuffer;
};

