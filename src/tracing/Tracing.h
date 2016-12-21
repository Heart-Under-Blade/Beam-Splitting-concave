#pragma once

#include "Beam.h"
#include "Particle.h"
#include "Intersection.h"

#include <float.h>
#include <vector>

#define MAX_BEAM_REFL_NUM 256

#define EPS_COS_90	1.7453292519943295769148298069306e-10	//cos(89.99999999)
#define EPS_COS_00	0.99999999998254670756866631966593		//1- cos(89.99999999)

// REF: поменять название на ~Splitter
class Tracing
{
public:
	Tracing(Particle *particle, const Point3f &startBeamDir, bool isOpticalPath,
			const Point3f &polarizationBasis, int interReflectionNumber);

	void RotateParticle(double beta, double gamma);

	virtual double BeamCrossSection(const Beam &/*beam*/) const { return 0.0; }

	virtual void SplitBeamByParticle(std::vector<Beam> &/*outBeams*/) {}

	virtual void SplitBeamByParticle(const std::vector<std::vector<int>> &tracks,
									 std::vector<Beam> &outBeams);

	double AreaOfBeam(const Beam &beam) const;

	double GetLightSurfaceArea() const;

protected:
	Particle *m_particle;			///< scattering particle (crystal)
	Facet *m_facets;
	Point3f m_polarizationBasis;	///<
	bool m_isOpticalPath;
	bool m_isArea;
	int m_interReflectionNumber;
	Beam m_initialBeam;				///< origin infinity beam

	Beam m_beamTree[MAX_BEAM_REFL_NUM];	///< tree of beams (works like stack)
	int m_treeSize;
	double m_lightSurfaceArea;

	const double FAR_ZONE_DISTANCE = 10000.0;
	const double LOW_ENERGY_LEVEL = 2e-12;

protected:

//	virtual void TraceInternalReflections(BeamInfo */*tree*/, int /*treesize*/,
//										  std::vector<Beam> &/*outBeams*/) {}

	void SetOpticalBeamParams_initial(int facetId, Beam &inBeam, Beam &outBeam);

	void SetBeam(Beam &beam, const Beam &other, const Point3f &dir, const Point3f &e,
				 const complex &coef1, const complex &coef2) const;

	bool Intersect(int facetIndex, const Beam& beam, Polygon &intersection) const;

	void SetPolygonByFacet(int facetId, Polygon &polygon) const;

	void CalcOpticalPathInternal(double cosIN, const Beam &incidentBeam, Beam &outBeam, Beam &inBeam) const;

	bool isTerminalBeam(const Beam &beam);

	void SplitExternalBeamByFacet(int facetId, Beam &inBeam, Beam &outBeam);

	void SplitInternalBeamByFacet(Beam &incidentBeam, int facetIndex,
								  Beam &inBeam, std::vector<Beam> &outBeams);

	void RotatePolarisationPlane(const Point3f &dir, const Point3f &facetNormal,
								 Beam &beam);

	void SetSloppingBeamParams_initial(const Point3f &beamDir, double cosIN, int facetId,
									   Beam &inBeam, Beam &outBeam);

	void SetNormalIncidenceBeamParams(double cosIN, const Beam &incidentBeam,
									  Beam &inBeam, Beam &outBeam);

	void SetSloppingIncidenceBeamParams(double cosIN, const Point3f &normal,
										Beam &incidentBeam, Beam &inBeam, Beam &outBeam,
										bool &isTrivialIncidence);

	void CalcLigthSurfaceArea(int facetId, const Beam &beam);

	void CalcOpticalPath_initial(Beam &inBeam, Beam &outBeam);

private:
	double ri_coef_re;
	double ri_coef_im;

private:
	double CalcNr(const double &cosIN) const;

	void SetTrivialIncidenceBeamParams(double cosIN, double Nr, const Point3f &normal, Point3f r0, double s,
									   const Beam &incidentBeam, Beam &inBeam, Beam &outBeam);

	void SetCompleteReflectionBeamParams(double cosIN, double Nr, const Beam &incidentBeam,
										 Beam &inBeam);

	void DivideBeamDirection(const Point3f &incidentDir, double cosIN, const Point3f &normal,
							 Point3f &reflDir, Point3f &refrDir) const;

	void SetOutputPolygon(__m128 *_output_points, int outputSize,
						  Polygon &polygon) const;

	bool ProjectToFacetPlane(const Polygon &polygon, const Point3f &dir,
							 const Point3f &normal, __m128 *_projection) const;

};
