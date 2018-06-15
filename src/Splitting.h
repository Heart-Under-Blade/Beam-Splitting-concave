#pragma once

#include "Beam.h"

class Splitting
{
public:
	Splitting(bool isOpticalPath);
	void ComputeRiParams(const complex &ri);

	void ComputeCosA(const Point3f &normal, const Point3f &incidentDir);
	void ComputeSplittingParams(const Point3f &dir, const Point3f &normal);

	bool IsCompleteReflection();
	bool IsNormalIncidence();
	bool IsIncident();

	double ComputeEffectiveReRi() const;

	double ComputeIncidentOpticalPath(const Point3f &direction,
									  const Point3f &facetPoint);
	double ComputeOutgoingOpticalPath(const Beam &beam);
	double ComputeSegmentOpticalPath(const Beam &beam,
									 const Point3f &facetPoint) const;

	void ComputeCRBeamParams(const Point3f &normal, const Beam &incidentBeam,
							 Beam &inBeam);

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

	complex m_ri;	//  refractive index
	double m_cRiRe;
	double m_cRiRe2;
	double m_cRiIm;

public:
	const double FAR_ZONE_DISTANCE = 10000.0; ///< distance from the center of coordinate system to the "far zone"

private:
	void ComputeCRJonesParams(complex &cv, complex &ch);

	void ComputeRegularJonesParams(const Point3f &normal,
								   const Beam &incidentBeam,
								   Beam &inBeam, Beam &outBeam);
	void ComputeInternalRefractiveDirection(const Vector3f &r,
											const Vector3f &normal, Vector3f &dir);

};
