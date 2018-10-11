#pragma once

#include "Tracks.h"
#include "Beam.h"
#include "Particle.h"
#include "Splitting.h"

#include <float.h>

//#define MAX_BEAM_REFL_NUM 32768
#define MAX_BEAM_REFL_NUM 65536
//#define MAX_BEAM_REFL_NUM 1048576

#define EPS_M_COS_90	-1.7453292519943295769148298069306e-10	//cos(89.99999999)

class Incidence;

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

struct OpticalPath
{
	double internal = 0;
	double external = 0;

	double GetTotal() { return internal + external; }
};

class Scattering
{
protected:
	Particle *m_particle;	///< scattering particle (crystal)
	int m_nActsMax;

	Splitting m_splitting;

	LightFacetChecker m_lightChecker;
	BeamFacetChecker m_beamChecker;

	Beam m_originBeam;
	Beam m_tracingBeams[MAX_BEAM_REFL_NUM];	///< tree of beams (works like stack)
	int m_treeSize;

	double m_incidentEnergy;

	const double EPS_BEAM_ENERGY = 2e-12;

public:
	Scattering(Particle *particle, Light *incidentLight, bool isOpticalPath,
			   int nActs);

	void ScatterLight(std::vector<Beam> &scaterredBeams);
	void RotateParticle(const Angle &angle);

	double GetIncedentEnergy() const;

	OpticalPath ComputeOpticalPath(const Beam &beam, const Point3f &startPoint,
								   std::vector<int> track);
	Particle *GetParticle() const;

protected:
	virtual void SplitLightToBeams(std::vector<Beam> &scatteredBeams) = 0;
	virtual void SplitBeams(std::vector<Beam> &scatteredBeams) = 0;

	virtual void FindVisibleFacets(const Beam &beam, FacetChecker &checker,
								   int begin, int end, Array<Facet*> &facets);

	bool ComputeOpticalBeamParams(Facet *facet, Beam beam);

	void SetPolygonByFacet(Facet *facet, Polygon &polygon) const;

	bool IsTerminalAct(const Beam &beam);

	Point3f ComputeBeamDirection(const Vector3f &oldDir, const Vector3f &normal,
								bool isIn1, bool isIn2);

	void ComputeFacetEnergy(const Vector3f &facetNormal,
							const Polygon &lightedPolygon);

	void PushBeamToTree(Beam &beam, Facet *facet, int level, bool isIn);

	void PushBeamToTree(Beam &beam, const Beam &oldBeam,
						const IdType &newId, Facet *facet, bool isIn);

	IdType RecomputeTrackId(const IdType &oldId, int facetId);
};
