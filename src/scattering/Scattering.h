#pragma once

#include "Tracks.h"
#include "Beam.h"
#include "Particle.h"
#include "Splitting.h"

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
 * @brief Produce a set of beams from a light that incident on a Particle.
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
	OpticalPath ComputeOpticalPath(const Beam &beam, const Point3f &startPoint,
								   std::vector<int> track);
protected:
	Particle *m_particle;	///< The crystal for a light scattering
	int m_maxActNo;			///< Maximal number of reflection/refraction acts

	Splitting m_splitting;

	LightFacetChecker m_lightChecker;
	BeamFacetChecker m_beamChecker;

	Beam m_originBeam;	///< The first beam that incident on a Particle.
						///< It has no boundaries (i.e. Polygon is empty)
	Beam m_propagatingBeams[MAX_BEAM_NUM];	///< Beams that are waiting for the next r/r act
	int m_treeSize;
	std::vector<Beam> *m_scatteredBeams;

	Array<Facet*> m_visibleFacets;

	double m_incidentEnergy;

	const double EPS_BEAM_ENERGY = 2e-12;

protected:
	virtual void SplitOriginBeam(std::vector<Beam> &scatteredBeams) = 0;

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

	void SplitSecondaryBeams(std::vector<Beam> &scatteredBeams);
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
};
