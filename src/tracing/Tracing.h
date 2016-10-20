#pragma once

#include "Beam.h"
#include "Particle.h"
#include "Intersection.h"

#include <float.h>
#include <vector>

#define MAX_BEAM_REFL_NUM 128

#define EPS_COS_89	1.7453292519943295769148298069306e-10	//cos(89.99999999)
#define EPS_COS_00	0.99999999998254670756866631966593		//1- cos(89.99999999)

struct BeamInfo
{
	Beam beam;
	int facetId;
	int dept;
};

class Tracing
{
public:
	Tracing(Particle *particle, const Point3f &startBeamDir, bool isOpticalPath,
			const Point3f &polarizationBasis, int interReflectionNumber);

	void RotateParticle(double beta, double gamma);

	virtual double BeamCrossSection(const Beam &/*beam*/) const {}

	virtual void SplitBeamByParticle(std::vector<Beam> &/*outBeams*/,
									 double &/*lightSurfaceSquare*/) {}

	virtual void SplitBeamByParticle(const std::vector<std::vector<int>> &tracks,
									 std::vector<Beam> &outBeams);

	double Square(const Beam &beam);
private:
	/**
	 * @brief The IncidenceCase enum
	 * Case of beam incidence to facet of particle
	 */
	enum class IncidenceCase : bool
	{
		Normal, Slopping
	};

protected:
	Particle *m_particle;			///< scattering particle (crystal)
	Point3f m_polarizationBasis;	///<
	bool m_isOpticalPath;
	int m_interReflectionNumber;
	Beam m_startBeam;

	const double FAR_ZONE_DISTANCE = 10000.0;
	const double LOW_ENERGY_LEVEL = 2e-12;

protected:

//	virtual void TraceInternalReflections(BeamInfo */*tree*/, int /*treesize*/,
//										  std::vector<Beam> &/*outBeams*/) {}

	void SetBeamsParamsExternal(int facetIndex, double cosIncident, Beam &inBeam, Beam &outBeam);


	void SetBeam(Beam &beam, const Beam &other, const Point3f &dir, const Point3f &e,
				 const complex &coef1, const complex &coef2) const;

	bool Intersect(int facetIndex, const Beam& originBeam, Beam &intersectBeam) const;

	void SplitBeamDirection(const Point3f &incidentDir, double cosIncident, const Point3f &normal,
								Point3f &refleDir, Point3f &refraDir) const;

	void SetBeamShapesByFacet(int facetIndex, Beam &inBeam, Beam &outBeam) const;

	void SetOutputBeam(__m128 *_output_points, int outputSize, Beam &outputBeam) const;

	bool ProjectToFacetPlane(const Beam& inputBeam, __m128 *_output_points,
							 __m128 _normal, int facetIndex) const;

	void CalcOpticalPathInternal(double Nr, const Beam &incidentBeam, Beam &outBeam, Beam &inBeam) const;

	void SplitBeam(const Beam &incidentBeam, Beam &inBeam, Beam &outBeam, double Nr,
				   IncidenceCase incidenceCase);

	void InvertBeamShapeOrder(Beam &outBeam, const Beam &inBeam);

	bool isEnough(const BeamInfo &info);

	void SplitExternalBeamByFacet(int facetIndex, double cosIncident,
								  Beam &inBeam, Beam &outBeam);

	void SplitInternalBeamByFacet(Beam &incidentBeam, int facetIndex,
								  Beam &inBeam, std::vector<Beam> &outBeams);
};
