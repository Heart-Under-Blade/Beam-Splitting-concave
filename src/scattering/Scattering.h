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
#define MAX_BEAM_NUM 65536
//#define MAX_BEAM_REFL_NUM 1048576

#define EPS_M_COS_90	-1.7453292519943295769148298069306e-10	//cos(89.99999999)
#define EPS_COS_90		1.7453292519943295769148298069306e-10	//cos(89.99999999)
#define EPS_COS_00		0.99999999998254670756866631966593		//1- cos(89.99999999)

class Incidence;

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
protected:
	Scattering(const complex &m_refractiveIndex,
			   int maxActNo, double minBeamEnergy = MIN_BEAM_NRG,
			   double farFresnelZone = FAR_FRESNEL_ZONE);

public:
	virtual ~Scattering();

	/**
	 * @brief Choose the Particle for scattering
	 * @param particle active particle
	 */
	void SetActiveParticle(Particle *particle);

	/**
	 * @brief Retrieve bunch of beams with different directions and other
	 * properties by scattering of the start beam by the Particle
	 * @param scaterredBeams result beams
	 */
	void Scatter(std::vector<Beam> *scatteredBeams);

	void Scatter(TrackNode */*trackTree*/,
				 std::vector<Beam> &/*scatteredBeams*/);

	virtual const Beam &ExtractShadowBeam();

	void SetStartBeam(const Vector3f &direction, const Vector3f &polarBasis);
	const Beam &GetStartBeam() const;

	double GetIncidentEnergy() const;

	/**
	 * @brief Computes optical path lenght of scattered beam
	 * @param beam scattered beam
	 * @param startPoint first point in optical path in shape of beam
	 * @param track set of indices of facets of the Particle
	 * @return internal and external optical path lenght
	 */
	OpticalPath ComputeOpticalPath(const BeamInfo &info,
								   const Point3f &startPoint);

	virtual void SelectOriginVisibleFacets(Array<Facet *> *facets);


	virtual double MeasureOpticalPath(const BeamInfo &info,
									  const Point3f sourcePoint);

	virtual double MeasureFullOpticalPath(const BeamInfo &info,
										  const Point3f sourcePoint);
//	double CrossSection(const Point3f &beamDir) const;

	void SetRefractiveIndex(const complex &ri);

	void SetFarFresnelZone(double value);
	double GetFarFresnelZone() const;

	complex m_refractiveIndex;	///< complex value of refractive index of the particle

protected:
	constexpr static const double MIN_BEAM_NRG = 2e-12;
	constexpr static const double FAR_FRESNEL_ZONE = 10000;

	TrackNode *m_trackTreeNode;
	bool m_hasTracks;

	int m_maxActNo;			///< Maximal number of reflection/refraction acts
	double m_minBeamEnergy;
	double m_farFresnelZone;

	Particle *m_particle;	///< The crystal for scattering
	Splitting *m_splitting;

	Incidence *m_incidence;
	RegularIncidence				m_regularIncidence;
	NormalIncidence					m_normalIncidence;
	TotalReflectionIncidence		m_totalReflectionIncidence;

	LightFacetChecker m_lightChecker;

	Beam m_shadowBeam;
	Beam m_startBeam; ///< The first beam that incident on a Particle. It has no boundaries (i.e. Polygon is empty)

	Stack_2p17<Beam> m_internalBeams; ///< Beams that are waiting for the next splitting (propagates only into the Particle)
	Stack_2p17<Beam>* m_secondaryBeams;
//	Beam m_activeBeams[MAX_ACTIVE_BEAM_NO];	///< Beams that are waiting for the next splitting
//	int m_nActiveBeams;
	std::vector<Beam> *m_scatteredBeams;

	Array<Facet*> m_workFacets;
	Array<Facet*> m_visibleFacets;

	double m_incidentEnergy;

	bool m_isBeamInside;

	void OrderVertices2f(std::vector<Point2f> &vertices,
						 Polygon &orderedPolygon);

	void RemoveDublicatedVertices2f(const std::vector<Point2f> &projected,
									std::vector<Point2f> &cleared);

protected:
	/**
	 * @brief Splits the original beam to external and internal beams
	 * @param externalBeams external beams
	 */
	virtual void SplitStartBeam() = 0;

	/**
	 * @brief Collect beams producted from splitting of the original beam and secondary beams
	 * @param scatteredBeams beams that leaved the Particle
	 */
	virtual void SplitSecondaryBeams() = 0;

	/**
	 * @brief Final handling of beam and throwing it out of the Particle
	 * @param beam throwed beam
	 * @param scatteredBeams output buffer
	 */
	virtual void ResolveExternalBeam(const Beam &beam) = 0;

	/**
	 * @brief Checks if beam has been to release out of Particle
	 * @param beam checked beam
	 * @return true if beam has been to release, false otherwise
	 */
	bool IsFinalAct(const Beam &beam);

	virtual bool isFinalFacet(int index, Array<Facet*> &facets);
	virtual void SelectVisibleFacets(const Beam &beam, Array<Facet*> *facets) = 0;
	virtual void ResolveBeams(Beam &parentBeam, Facet *facet);

	void Reset();

	void ComputeSplittingParams(const Vector3f &dir, const Vector3f &normal,
								bool isInside);

	void SplitBeamByVisibleFacets(Beam &beam);

	Array<Facet*> *FindVisibleFacets(const Beam &beam, FacetChecker &checker,
									 Array<Facet*> &facets,
									 Array<Facet*> *visibleFacets);

	void ComputeOpticalBeamParams(const Facet *facet, Beam beam,
								  const Polygon &resultShape);

	void SetPolygonByFacet(Facet *facet, Polygon &polygon) const;

	Point3f ComputeBeamDirection(const Vector3f &oldDir, const Vector3f &normal,
								 bool isIn1, bool isIn2);

	void ComputeFacetEnergy(const Vector3f &facetNormal,
							const Polygon &lightedPolygon);

	void StackBeam(const Beam &beam);
	void SetIncidence();

private:
	void SetShadowBeamsParams();
};

//int Scattering::MAX_ACTIVE_BEAM_NUM ;
