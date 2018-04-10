#pragma once

#include "Beam.h"
#include "Particle.h"
#include "Intersection.h"

#include <float.h>

//#define MAX_BEAM_REFL_NUM 32768
#define MAX_BEAM_REFL_NUM 65536
//#define MAX_BEAM_REFL_NUM 1048576

#define EPS_M_COS_90	-1.7453292519943295769148298069306e-10	//cos(89.99999999)
#define EPS_COS_90		1.7453292519943295769148298069306e-10	//cos(89.99999999)
#define EPS_COS_00		0.99999999998254670756866631966593		//1- cos(89.99999999)

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

	Light *m_incidentLight;
	Point3f m_incidentDir;
	Point3f m_polarBasis;

	bool m_isOpticalPath;
	int m_nActs;

//	std::vector<Beam> m_beamTree;
	Beam m_beamTree[MAX_BEAM_REFL_NUM];	///< tree of beams (works like stack)
	int m_treeSize;
	double m_incommingEnergy;

	const double FAR_ZONE_DISTANCE = 10000.0;
	const double LOW_ENERGY_LEVEL = 2e-12;

private:
	double m_cRiRe;
	double m_cRiRe2;
	double m_cRiIm;
	complex m_refrIndex;

public:
	Scattering(Particle *particle, Light *incidentLight, bool isOpticalPath,
			   int nActs);

	virtual void ScatterLight(double /*beta*/, double /*gamma*/, std::vector<Beam> &/*scaterredBeams*/) {}
	virtual void ScatterLight(double beta, double gamma, const std::vector<std::vector<int>> &tracks,
									 std::vector<Beam> &scaterredBeams);

	double GetIncomingEnergy() const;

	double ComputeInternalOpticalPath(const Beam &beam, const std::vector<int> &track);
//	double CrossSection(const Point3f &beamDir) const;

protected:
	void SetBeamOpticalParams(int facetId, Beam &inBeam, Beam &outBeam);

	void SetBeam(Beam &beam, const Beam &other, const Point3f &dir,
				 const complex &coef1, const complex &coef2) const;

	void Difference(const Polygon &subject, const Point3f &subjNormal,
					const Polygon &clip, const Point3f &clipNormal,
					const Point3f &clipDir, Polygon *difference, int &resultSize) const;

	bool Intersect(int facetId, const Beam& beam, Polygon &intersection) const;

	void SetPolygonByFacet(int facetId, Polygon &polygon) const;

	void CalcOpticalPath(double cosIN, const Beam &incidentBeam, Beam &inBeam,
						 Beam &outBeam) const;

	bool IsTerminalAct(const Beam &beam);

	void SplitLightToBeams(int facetId, Beam &inBeam, Beam &outBeam);

	void SplitSecondaryBeams(Beam &incidentBeam, int facetID,
								  Beam &inBeam, std::vector<Beam> &outBeams);

	void RotatePolarisationPlane(const Point3f &dir, const Point3f &facetNormal,
								 Beam &beam);

	void SetRegularBeamParamsExternal(const Point3f &beamDir, double cosIN, int facetId,
									   Beam &inBeam, Beam &outBeam);

	void SetNormalIncidenceBeamParams(double cosIN, const Beam &incidentBeam,
									  Beam &inBeam, Beam &outBeam);

	void SetRegularIncidenceBeamParams(double cosIN, const Point3f &normal,
									   Beam &incidentBeam,
									   Beam &inBeam, Beam &outBeam,
									   bool &isTrivialIncidence);

	void CalcFacetEnergy(int facetId, const Polygon &lightedPolygon);

	void CalcOpticalPathForLight(Beam &inBeam, Beam &outBeam);

	void PushBeamToTree(Beam &beam, int facetId, int level, Location location);
	void PushBeamToTree(Beam &beam, int facetId, int level);
	void PushBeamToTree(Beam &beam);

	void SetBeamID(Beam &beam);

private:
	double CalcReRI(const double &cosIN) const;

	void SetRegularBeamParams(double cosIN, double reRI, const Point3f &normal,
									   Point3f r0, double s, const Beam &incidentBeam,
									   Beam &inBeam, Beam &outBeam);

	void SetCRBeamParams(double cosIN, double Nr, const Beam &incidentBeam,
										 Beam &inBeam);

	void SplitDirection(const Point3f &incidentDir, double cosIN, const Point3f &normal,
							 Point3f &reflDir, Point3f &refrDir) const;

	Point3f ChangeBeamDirection(const Vector3f &oldDir, const Vector3f &normal,
								Location loc);

	void SetOutputPolygon(__m128 *_output_points, int outputSize,
						  Polygon &polygon) const;

	bool ProjectToFacetPlane(const Polygon &polygon, const Point3f &dir,
							 const Point3f &normal, __m128 *_projection) const;
};
