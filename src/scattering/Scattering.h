#pragma once

#include "Tracks.h"
#include "Beam.h"
#include "Particle.h"
#include "Splitting.h"

#include "CompleteReflectionIncidence.h"
#include "NormalIncidence.h"

#include <float.h>

//#define MAX_BEAM_REFL_NUM 32768
#define MAX_BEAM_NUM 65536
//#define MAX_BEAM_REFL_NUM 1048576

#define EPS_M_COS_90	-1.7453292519943295769148298069306e-10	//cos(89.99999999)

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
 * @brief Transform light that incidents on a Particle to set of beams.
 */
class Scattering
{
public:
	Scattering(Particle *particle, const Light &incidentLight, int maxActNo);
	virtual ~Scattering();

	/**
	 * @brief Tranform incident light into beams after throwing the Particle
	 * @param scaterredBeams result beams
	 */
	void ScatterLight(std::vector<Beam> &scaterredBeams);

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
protected:
	Particle *m_particle;	///< The crystal for a light scattering
	int m_maxActNo;			///< Maximal number of reflection/refraction acts

	Splitting m_splitting;

	Incidence *m_incidence;
	RegularIncidence				m_regularIncidence;
	NormalIncidence					m_normalIncidence;
	CompleteReflectionIncidence		m_completeReflectionIncidence;

	LightFacetChecker m_lightChecker;
	BeamFacetChecker m_beamChecker;

	Beam m_originalBeam;	///< The first beam that incident on a Particle.
							///< It has no boundaries (i.e. Polygon is empty)
	Beam m_propagatingBeams[MAX_BEAM_NUM];	///< Beams that are waiting for the next r/r act
	int m_treeSize;
	std::vector<Beam> *m_scatteredBeams;

	Array<Facet*> m_visibleFacets;

	double m_incidentEnergy;
	const double EPS_BEAM_ENERGY = 2e-12;

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
	virtual bool IsTerminalAct(const Beam &beam);
	virtual bool isTerminalFacet(int index, Array<Facet*> &facets);
	virtual void SelectVisibleFacets(const Beam &beam, Array<Facet*> &facets) = 0;
	virtual void PushBeamsToBuffer(Beam &parentBeam, Facet *facet,
								   bool hasOutBeam);

	void ComputeSplittingParams(const Vector3f &dir, const Vector3f &normal,
								bool isInside);

	void SplitBeamByVisibleFacets(Beam &beam);

	void FindVisibleFacets(const Beam &beam, FacetChecker &checker,
						   int begin, int end, Array<Facet*> &facets);

	bool ComputeOpticalBeamParams(Facet *facet, Beam beam);

	void SetPolygonByFacet(Facet *facet, Polygon &polygon) const;

	Point3f ComputeBeamDirection(const Vector3f &oldDir, const Vector3f &normal,
								 bool isIn1, bool isIn2);

	void ComputeFacetEnergy(const Vector3f &facetNormal,
							const Polygon &lightedPolygon);

	void PushBeamToTree(Beam &beam);
	void CreateOriginBeam(const Vector3f &dir, const Vector3f &basis);
	void SetIncidence();
};
