#pragma once

#include "Beam.h"

class Splitting
{
public:
	Splitting(bool isOpticalPath, Light *incidentLight);
	void ComputeRiParams(const complex &ri);

	void ComputeSplittingParams(const Point3f &dir, const Point3f &normal);

	bool IsCompleteReflection();
	bool IsNormalIncidence();
	bool IsIncident();

	double ComputeEffectiveReRi() const;

	double ComputeIncidentOpticalPath(const Point3f &facetPoint);
	double ComputeScatteredOpticalPath(const Beam &beam);
	double ComputeSegmentOpticalPath(const Beam &beam,
									 const Point3f &facetPoint) const;

	void ComputeCRBeamParams(const Point3f &normal, const Beam &incidentBeam,
							 Beam &inBeam);

	void ComputeCosA(const Point3f &normal, const Point3f &incidentDir);


	void ComputeNormalBeamParams(const Beam &incidentBeam,
								 Beam &inBeam, Beam &outBeam);
	void ComputeNormalBeamParamsExternal(const Light &incidentLight,
										 Beam &inBeam, Beam &outBeam);

	void ComputeRegularBeamsParams(const Point3f &normal,
								   const Beam &incidentBeam,
								   Beam &inBeam, Beam &outBeam);
	void ComputeRegularBeamParamsExternal(const Point3f &facetNormal,
										  Beam &incidentBeam,
										  Beam &inBeam, Beam &outBeam);

	Point3f ChangeBeamDirection(const Vector3f &oldDir, const Vector3f &normal,
								Location oldLoc, Location loc);
private:
	Point3f r;
	double reRiEff;
	double s;
	double cosA;
	bool m_isOpticalPath;
	Light *m_incidentLight;

	const double FAR_ZONE_DISTANCE = 10000.0;

	complex m_ri;	//  refractive index
	double m_cRiRe;
	double m_cRiRe2;
	double m_cRiIm;

private:
	void ComputeCRJonesParams(complex &cv, complex &ch);

	void ComputeRegularJonesParams(const Point3f &normal,
								   const Beam &incidentBeam,
								   Beam &inBeam, Beam &outBeam);
	void ComputeInternalRefractiveDirection(const Vector3f &r,
											const Vector3f &normal, Vector3f &dir);

};
