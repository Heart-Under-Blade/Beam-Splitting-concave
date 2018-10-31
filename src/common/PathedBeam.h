#pragma once

#include "Beam.h"

/**
 * @brief A beam with an optical path
 */
class PathedBeam : public Beam
{
public:
	PathedBeam();

	void Clear();
	void AddOpticalPath(double path);
	double ComputeOutgoingOpticalPath() const;
	double ComputeIncidentOpticalPath(const Point3f &direction,
									  const Point3f &facetPoint) const;

	double ComputeSegmentOpticalPath(const double &reRiEff,
									 const Point3f &facetPoint) const;

	PathedBeam & operator = (const PathedBeam &other);
	PathedBeam & operator = (PathedBeam &&other);

public:
	// REF: перенести в private
	double opticalPath;	///< Optical path of beam
	double front;		///< Current position of phase front from Ax+By+Cz+D=0 (where D is front)

#ifdef _DEBUG // DEB
	std::vector<double> ops;
#endif

private:
	const double FAR_ZONE_DISTANCE = 10000.0; ///< Distance from the center of coordinate system to the "far zone"
};
