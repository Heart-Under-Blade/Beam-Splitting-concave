#pragma once

#include "Scattering.h"

/** NOTE: пучки выходят со случайно ориентированным порядком вершин */
class ScatteringNonConvex : public Scattering
{
public:
	ScatteringNonConvex(Particle *particle, Light *incidentLight,
						bool isOpticalPath, int nActs);
private:
	void SortFacets_faster(const Point3f &beamDir, Array<Facet*> &facets) const;
	int FindClosestVertex(const Polygon &facet, const Point3f &beamDir) const;
	void CutBeamByFacet(Facet *facet, Beam &beam,
						PolygonArray &result);

	double CalcMinDistanceToFacet(Polygon *facet, const Point3f &beamDir);
	void SortFacets(const Point3f &beamDir, Array<Facet*> &facets); ///< use 'Fast sort' algorithm

	void CutExternalBeam(const Beam &beam, std::vector<Beam> &scaterredBeams);

	void FindVisibleFacetsForBeam(const Beam &beam, Array<Facet *> &facets);

	void SelectVisibleFacets(const Beam &beam, Array<Facet *> &facets);
	void SelectVisibleFacetsForLight(Array<Facet*> &facets);

	bool IntersectLightWithFacet(const Array<Facet*> &facets, int nCheckedFacets,
								 PolygonArray &resFacets);


	bool IsOutgoingBeam(Beam &incidentBeam);

	int FindFacetId(int facetId, const Array<Facet*> &arr);

	void TraceFirstBeamFixedFacet(int facetID, bool &isIncident);

	void PushBeamsToTree(Facet *facet, Splitting &splitting,
						 const PolygonArray &polygons);

	void SplitBeamByVisibleFacets(Beam &beam);

	void PushBeamsToBuffer(Facet *facet, const Beam &beam, bool hasOutBeam,
						   Beam &inBeam, Beam &outBeam, std::vector<Beam> &passed);

	void CutPolygonByFacets(const Polygon &pol,
							const Array<Facet*> &facets, int size,
							const Vector3f &polNormal, const Vector3f &clipNormal,
							const Vector3f &dir,
							PolygonArray &pols);

	void PushBeamPartsToTree(const Beam &beam,
							 const PolygonArray &parts);
protected:
	void SplitBeams(std::vector<Beam> &scaterredBeams) override;
	void SplitLightToBeams(std::vector<Beam> &/*scatteredBeams*/) override;

	bool FindRestOfBeamShape(Facet *facet, Beam &beam);

private:
	PolygonArray m_polygonBuffer;
};

