#pragma once

#include "Scattering.h"

/** NOTE: пучки выходят со случайно ориентированным порядком вершин */

/**
 * @brief Provide methods for non-convex type of particles.
 */
class ScatteringNonConvex : public Scattering
{
public:
	ScatteringNonConvex(Particle *particle, const Light &incidentLight, int maxActNo);
	~ScatteringNonConvex();

protected:
	void SplitOriginalBeam(std::vector<Beam> &externalBeams) override;

	void ReleaseBeam(Beam &beam) override;
	bool IsTerminalAct(const Beam &beam) override;
	bool isTerminalFacet(int index, Array<Facet*> &facets) override;
	void PushBeamsToBuffer(Beam &parentBeam, Facet *facet,
						   bool hasOutBeam) override;
	void SelectVisibleFacets(const Beam &beam, Array<Facet*> &facets) override;

	void PushBeamToBuffer(Beam &beam, const PolygonArray &beamParts,
						  std::vector<Beam> &scatteredBeams);

private:
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

	bool FindRestOfBeamShape(Facet *facet, const Beam &beam, PolygonArray &rest);

	double CalcMinDistanceToFacet(Polygon *facet, const Point3f &beamDir);
	void SortFacets(const Point3f &beamDir, Array<Facet*> &facets); ///< use 'Fast sort' algorithm

	bool FindLightedFacetPolygon(const Array<Facet*> &facets, int nCheckedFacets,
								 PolygonArray &pols);

	void PushBeamsToTree(Facet *facet, BeamPair<Beam> &beams,
						 const PolygonArray &polygons,
						 std::vector<Beam> &scatteredBeams);

	void CutPolygonByFacets(const Polygon &pol,
							const Array<Facet*> &facets, int size,
							const Vector3f &polNormal, const Vector3f &clipNormal,
							const Vector3f &dir, PolygonArray &pols);

	void PushBeamPartsToBuffer(const Beam &beam, const PolygonArray &parts);
private:
	bool m_isDivided;
	PolygonArray m_intersectionBuffer;	///< Buffer for result of Polygon intersection functions
	PolygonArray m_differenceBuffer;	///< Buffer for result of Polygon differencefunctions
};

