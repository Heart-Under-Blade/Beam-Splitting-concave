#pragma once

#include "Tracks.h"
#include "Beam.h"
#include "Particle.h"
#include "Intersection.h"
#include "Splitting.h"
#include "TrackTree.h"

#include "TotalReflectionIncidence.h"
#include "NormalIncidence.h"

#include "TotalReflectionIncidence.h"
#include "NormalIncidence.h"

#include <float.h>

//#define MAX_BEAM_REFL_NUM 32768
#define MAX_BEAM_REFL_NUM 65536
//#define MAX_BEAM_REFL_NUM 1048576

#define EPS_M_COS_90	-1.7453292519943295769148298069306e-10	//cos(89.99999999)

/**
 * @brief The BeamTree struct
 * Tree of beams (works like stack).
 */
struct BeamTree
{
	Beam tree[MAX_BEAM_REFL_NUM];
	int size = 0;

//	void Push(const Beam &beam)
//	{
//		tree[size++] = beam;
//	}
};

class Scattering
{
public:
	// REF: убрать в протектед потом
	Particle *m_particle;			///< scattering particle (crystal)

protected:
	Facet *m_facets;
	Splitting m_splitting;
	Light *m_incidentLight;

	Vector3f m_incidentDir;
	Vector3f m_polarBasis;
	int m_nActs;

	Beam m_beamTree[MAX_BEAM_REFL_NUM];	///< tree of beams (works like stack)
	int m_treeSize;
	double m_incidentEnergy;

	const double EPS_BEAM_ENERGY = 2e-12;

public:
	Scattering(Particle *particle, const Light &incidentLight, int maxActNo,
			   const complex &m_refractiveIndex);
	virtual ~Scattering();

	virtual void ScatterLight(std::vector<Beam> &/*scaterredBeams*/) {}
	virtual void ScatterLight(const std::vector<std::vector<int>> &tracks,
									 std::vector<Beam> &scaterredBeams);
	void ScatterLight(TrackNode */*trackTree*/,
					  std::vector<Beam> &/*scatteredBeams*/);

	virtual void ExtractShadowBeam(std::vector<Beam> &scaterredBeams);

	double GetIncedentEnergy() const;

	virtual double MeasureOpticalPath(const Beam &beam, const Point3f sourcePoint,
									 const std::vector<int> &track);

	virtual double MeasureFullOpticalPath(const Beam &beam, const Point3f sourcePoint,
										  const std::vector<int> &track);
//	double CrossSection(const Point3f &beamDir) const;

	void SetSplitting(Particle *p);

protected:
	void SetIncidentBeamOpticalParams(unsigned facetId, Beam &inBeam, Beam &outBeam);

	void Difference(const Polygon &subject, const Vector3f &subjNormal,
					const Polygon &clip, const Vector3f &clipNormal,
					const Vector3f &clipDir, PolygonArray &difference) const;

	void Intersect(int facetId, const Beam& beam, Polygon &intersection) const;

	void SetPolygonByFacet(int facetId, Polygon &polygon) const;

	bool IsTerminalAct(const Beam &beam);

	virtual const Beam &SetShadowBeam();

	double GetIncidentEnergy() const;

	void SplitLightToBeams(int facetId, Beam &inBeam, Beam &outBeam);

	void ComputePolarisationParams(const Vector3f &dir,
								   const Vector3f &facetNormal, Beam &beam);

	void ComputeFacetEnergy(int facetId, const Polygon &lightedPolygon);
	complex m_refractiveIndex;	///< complex value of refractive index of the particle

protected:
	int m_maxActNo;			///< Maximal number of reflection/refraction acts

	void PushBeamToTree(Beam &beam, int facetId, int level, Location location);

	Incidence *m_incidence;
	RegularIncidence				m_regularIncidence;
	NormalIncidence					m_normalIncidence;
	TotalReflectionIncidence		m_totalReflectionIncidence;


	IdType RecomputeTrackId(const IdType &oldId, int facetId);
	Beam m_shadowBeam;
	Beam m_originalBeam;	///< The first beam that incident on a Particle.
							///< It has no boundaries (i.e. Polygon is empty)
	Beam m_propagatingBeams[MAX_BEAM_REFL_NUM];	///< Beams that are waiting for the next r/r act
	std::vector<Beam> *m_scatteredBeams;

	Array<Facet*> m_workFacets;
	Array<Facet*> m_visibleFacets;

	void OrderVertices2f(std::vector<Point2f> &vertices,
						 Polygon &orderedPolygon);

	void ProjectParticleToXY(std::vector<Point2f> &projected);

	void RemoveDublicatedVertices2f(const std::vector<Point2f> &projected,
								  std::vector<Point2f> &cleared);

private:
	void SetOutputPolygon(__m128 *_output_points, int outputSize,
						  Polygon &polygon) const;

	bool ProjectToFacetPlane(const Polygon &polygon, const Vector3f &dir,
							 const Point3f &normal, __m128 *_projection) const;

protected:
	TrackNode *m_trackTreeNode;
	bool m_hasTracks;

protected:
	/**
	 * @brief Splits the original beam to external and internal beams
	 * @param externalBeams external beams
	 */
	virtual void SplitOriginalBeam(std::vector<Beam> &externalBeams) = 0;

	/**
	 * @brief Collect beams producted from splitting of the original beam and secondary beams
	 * @param scatteredBeams beams that leaved the Particle
	 */
	void SplitSecondaryBeams(std::vector<Beam> &scatteredBeams);

	/**
	 * @brief Final handling of beam and throwing it out of the Particle
	 * @param beam throwed beam
	 * @param scatteredBeams output buffer
	 */
	virtual void ReleaseBeam(Beam &beam);

	/**
	 * @brief Checks if beam has been to release out of Particle
	 * @param beam checked beam
	 * @return true if beam has been to release, false otherwise
	 */
	virtual bool IsFinalAct(const Beam &beam);
	virtual bool isFinalFacet(int index, Array<Facet*> &facets);
	virtual void SelectVisibleFacets(const Beam &beam, Array<Facet*> &facets) = 0;
	virtual void PushBeamsToBuffer(Beam &parentBeam, Facet *facet,
								   bool hasOutBeam);

	void ComputeSplittingParams(const Vector3f &dir, const Vector3f &normal,
								bool isInside);

	void SplitBeamByVisibleFacets(Beam &beam);

	void FindVisibleFacets(const Beam &beam, FacetChecker &checker,
						   Array<Facet *> &facets,
						   Array<Facet*> &visibleFacets);

	bool ComputeOpticalBeamParams(Facet *facet, Beam beam,
								  const Polygon &resultShape);

	void SetPolygonByFacet(Facet *facet, Polygon &polygon) const;

	Point3f ComputeBeamDirection(const Vector3f &oldDir, const Vector3f &normal,
								 bool isIn1, bool isIn2);

	void ComputeFacetEnergy(const Vector3f &facetNormal,
							const Polygon &lightedPolygon);

	void PushBeamToTree(Beam &beam);
	void CreateOriginalBeam(const Vector3f &dir, const Vector3f &basis);
	void SetIncidence();
};
