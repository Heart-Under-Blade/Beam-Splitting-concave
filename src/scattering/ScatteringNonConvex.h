#pragma once

#include "Scattering.h"

/** NOTE: пучки выходят со случайно ориентированным порядком вершин */

/**
 * @brief Provide methods for non-convex type of particles.
 */
class ScatteringNonConvex : public Scattering
{
public:
	ScatteringNonConvex(const complex &refractiveIndex, int maxActNo,
						double minBeamEnergy = MIN_BEAM_NRG,
						double farFresnelZone = FAR_FRESNEL_ZONE);
	~ScatteringNonConvex();

	void Scatter(TrackNode *trackTree, std::vector<Beam> &scatteredBeams);

	double MeasureOpticalPath(const BeamInfo &info,
							  const Point3f sourcePoint) override;

protected:
	void SplitStartBeam() override;
	void SplitSecondaryBeams() override;
	void ResolveExternalBeam(const Beam &beam) override;
	bool isFinalFacet(int index, Array<Facet*> &facets) override;
	void ResolveBeams(Beam &parentBeam, Facet *facet) override;
	void SelectVisibleFacets(const Beam &beam, Array<Facet*> *visibleFacets) override;

private:
	/**
	 * @brief Sorts facets of the Particle by distance from the beginnig of
	 * latest beam splitting in descending order.
	 * @param beamDir beam direction
	 * @param facets unordered facets
	 */
	void SortFacetsByDistance(const Vector3f &beamDir, Array<Facet *> *facets) const;

	/**
	 * @brief Finds vertex of the Facet polygon that located more near than
	 * others to beginnig of the Beam.
	 * @param facet facet polygon
	 * @param beamDir beam direction
	 * @return index of the closest vertex in the beam polygon
	 */
	int FindClosestVertex(const Polygon &facet, const Point3f &beamDir) const;

	PolygonResult FindRestOfBeamShape(Facet *facet, const Beam &beam,
									  PolygonStack &rest);

	double CalcMinDistanceToFacet(Polygon *facet, const Point3f &beamDir);
	void SortFacets(const Point3f &beamDir, Array<Facet*> &facets); ///< use 'Fast sort' algorithm

	bool FindVisibleFacetPart(const Array<Facet*> &facets, int facetIndex,
								 PolygonStack &pols);

	void ResolveBeamParts(Facet *facet, const PolygonStack &parts);

	void CutPolygonByFacets(const Polygon &pol,
							const Array<Facet*> &facets, int size,
							const Vector3f &polNormal, const Vector3f &clipNormal,
							const Vector3f &dir, PolygonStack &pols);

	void StackBeamParts(Beam beam, const PolygonStack &parts);
	void SplitExternalBeam(const Beam &beam);
	void ReleaseFinalBeams(Beam beam);

private:
	bool m_isDivided;
	PolygonStack m_splittedBeams;	///< Buffer for result of Polygon intersection functions
	PolygonStack m_beamParts;	///< Buffer for result of Polygon differencefunctions
	BeamFacetChecker m_beamChecker;
	Stack_2p17<Beam> m_externalBeams;
};



