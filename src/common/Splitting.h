#pragma once

#include "Beam.h"

class Splitting
{
public:
	Splitting(bool isOpticalPath);
	void ComputeRiParams(const complex &ri);

	void ComputeCosA(const Point3f &normal, const Point3f &incidentDir);
	void ComputeSplittingParams(const Point3f &dir, const Point3f &normal,
								bool isInside);

	void ComputeReflectedDirection(Vector3f &dir) const
	{
		dir = r - m_normal;
		Point3f::Normalize(dir); // REF, OPT: нужно ли это нормализовать всегда?
	}

	void ComputeRefractedDirection(Vector3f &dir) const
	{
		dir = r/sqrt(s) + m_normal;
		Point3f::Normalize(dir); // REF, OPT: нужно ли это нормализовать всегда?
	}

	bool IsCompleteReflection();
	bool IsNormalIncidence();
	bool IsIncident();

	void SetBeams(const Polygon &beamShape);
	void SetNormal(const Point3f &normal);

	double ComputeEffectiveReRi() const;

	double ComputeIncidentOpticalPath(const Point3f &direction,
									  const Point3f &facetPoint);
	double ComputeOutgoingOpticalPath(const Beam &beam);
	double ComputeSegmentOpticalPath(const Beam &beam,
									 const Point3f &facetPoint) const;

	complex GetRi() const;

public:
	Point3f r;
	double reRiEff;
	double s;
	bool m_isOpticalPath;

	complex m_ri;	//  refractive index
	double m_cRiRe;
	double m_cRiRe2;
	double m_cRiIm;

	Beam inBeam;
	Beam outBeam;

	Point3f m_normal;

public:
	double cosA;
	const double FAR_ZONE_DISTANCE = 10000.0; ///< distance from the center of coordinate system to the "far zone"
};
