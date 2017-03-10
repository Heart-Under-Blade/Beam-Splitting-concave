#pragma once

#include "Tracing.h"

class TracingConvex : public Tracing
{
public:
	TracingConvex(Particle *particle, const Point3f &incidentBeamDir,
				  bool isOpticalPath, const Point3f &polarizationBasis,
				  int interReflectionNumber);

	void SplitBeamByParticle(double beta, double gamma, std::vector<Beam> &outBeams) override;
	void SplitBeamByParticle(double, double, const std::vector<std::vector<int>> &,
							 std::vector<Beam> &) override; ///> for predefined trajectories

protected:
	void TraceInternalReflections(std::vector<Beam> &outBeams);
};
