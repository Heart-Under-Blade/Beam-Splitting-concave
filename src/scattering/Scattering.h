#pragma once

#include "Tracks.h"
#include "Beam.h"
#include "Particle.h"
#include "Intersection.h"
#include "Splitting.h"

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
	Splitting m_splitting;
	Light *m_incidentLight;

	Vector3f m_incidentDir;
	Vector3f m_polarBasis;
	int m_nActs;

	Beam m_propagatingBeams[MAX_BEAM_REFL_NUM];	///< tree of beams (works like stack)
	int m_treeSize;
	double m_incidentEnergy;

	const double EPS_BEAM_ENERGY = 2e-12;

public:
	Scattering(Particle *particle, Light *incidentLight, bool isOpticalPath,
			   int nActs);

	virtual void ScatterLight(std::vector<Beam> &scaterredBeams) = 0;
	virtual void ScatterLight(const std::vector<std::vector<int>> &tracks,
							  std::vector<Beam> &scaterredBeams) = 0;
	void RotateParticle(const Angle &angle);

	double GetIncedentEnergy() const;

	OpticalPath ComputeOpticalPath(const Beam &beam, const Point3f &startPoint,
								   std::vector<int> track);
//	double CrossSection(const Point3f &beamDir) const;

	Particle *GetParticle() const;

protected:
	void SetIncidentBeamOpticalParams(Facet *facet, Beam &inBeam, Beam &outBeam);

	void Difference(const Polygon &subject, const Vector3f &subjNormal,
					const Polygon &clip, const Vector3f &clipNormal,
					const Vector3f &clipDir, PolygonArray &difference) const;

	bool IncidentBeamToFacet(Facet *facet, const Beam& beam,
							 Polygon &intersection) const;

	void SetPolygonByFacet(Facet *facet, Polygon &polygon) const;

	bool IsTerminalAct(const Beam &beam);

	void SplitLightToBeams(Facet *facet, Beam &inBeam, Beam &outBeam);

	void ComputePolarisationParams(const Vector3f &dir,
								   const Vector3f &facetNormal, Beam &beam);

	void ComputeFacetEnergy(const Vector3f &facetNormal,
							const Polygon &lightedPolygon);


	void PushBeamToTree(Beam &beam, Facet *facet, int level, Location location);


	IdType RecomputeTrackId(const IdType &oldId, int facetId);

private:
	void SetOutputPolygon(__m128 *_output_points, int outputSize,
						  Polygon &polygon) const;

	bool ProjectToFacetPlane(const Polygon &polygon, const Vector3f &dir,
							 const Point3f &normal, __m128 *_projection) const;

};
