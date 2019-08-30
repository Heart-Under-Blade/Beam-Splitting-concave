#pragma once

#include "Tracks.h"
#include "Beam.h"
#include "Particle.h"
#include "Splitting.h"
#include "TrackTree.h"

#include "TotalReflectionIncidence.h"
#include "NormalIncidence.h"

#include "TotalReflectionIncidence.h"
#include "NormalIncidence.h"

#include <float.h>

//#define MAX_BEAM_REFL_NUM 32768
#define MAX_BEAM_NUM 65536
//#define MAX_BEAM_REFL_NUM 1048576

#define EPS_M_COS_90	-1.7453292519943295769148298069306e-10	//cos(89.99999999)
#define EPS_COS_90		1.7453292519943295769148298069306e-10	//cos(89.99999999)
#define EPS_COS_00		0.99999999998254670756866631966593		//1- cos(89.99999999)

class Incidence;

/**
 * @brief Buffer for beams
 */
class BeamStack
{
	Beam arr[MAX_BEAM_NUM];
	int size = 0;

	void Push(const Beam &beam)
	{
		arr[size++] = beam;
	}
};

struct OpticalPath
{
	double internal = 0;
	double external = 0;

	double GetTotal() { return internal + external; }
};

/**
 * @brief Transform a light that incidents on a Particle to a set of beams.
 */
class Scattering
{
public:
	Scattering(Particle *particle, const Light &incidentLight, int maxActNo,
			   const complex &m_refractiveIndex);
	virtual ~Scattering();

	/**
	 * @brief Tranform incident light into beams after throwing the Particle
	 * @param scaterredBeams result beams
	 */
	void ScatterLight(std::vector<Beam> &scatteredBeams);

	void ScatterLight(TrackNode */*trackTree*/,
					  std::vector<Beam> &/*scatteredBeams*/);

	virtual const Beam &SetShadowBeam();

	double GetIncidentEnergy() const;

	/**
	 * @brief Computes optical path lenght of scattered beam
	 * @param beam scattered beam
	 * @param startPoint first point in optical path in shape of beam
	 * @param track set of indices of facets of the Particle
	 * @return internal and external optical path lenght
	 */
	OpticalPath ComputeOpticalPath(const Beam &beam, const Point3f &startPoint,
								   std::vector<int> track);

	virtual void SelectOriginVisibleFacets(Array<Facet*> &facets);

	complex m_refractiveIndex;	///< complex value of refractive index of the particle

protected:
	int m_maxActNo;			///< Maximal number of reflection/refraction acts
	Particle *m_particle;	///< The crystal for a light scattering

	Splitting m_splitting;

	Incidence *m_incidence;
	RegularIncidence				m_regularIncidence;
	NormalIncidence					m_normalIncidence;
	TotalReflectionIncidence		m_totalReflectionIncidence;

	LightFacetChecker m_lightChecker;
	BeamFacetChecker m_beamChecker;

	Beam m_shadowBeam;
	Beam m_originalBeam;	///< The first beam that incident on a Particle.
							///< It has no boundaries (i.e. Polygon is empty)
	Beam m_propagatingBeams[MAX_BEAM_NUM];	///< Beams that are waiting for the next r/r act
	int m_treeSize;
	std::vector<Beam> *m_scatteredBeams;

	Array<Facet*> m_workFacets;
	Array<Facet*> m_visibleFacets;

	double m_incidentEnergy;
	const double EPS_BEAM_ENERGY = 2e-12;

	void OrderVertices2f(std::vector<Point2f> &vertices,
						 Polygon &orderedPolygon);

	void ProjectParticleToXY(std::vector<Point2f> &projected);

	void RemoveDublicatedVertices2f(const std::vector<Point2f> &projected,
								  std::vector<Point2f> &cleared);

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
