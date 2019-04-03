#pragma once

#include "Tracing.h"

class TracingConvex : public Tracing
{
public:
	TracingConvex(Particle *particle, Light *incidentLight,
				  bool isOpticalPath, int interReflectionNumber);

	void SplitBeamByParticle(double beta, double gamma, std::vector<Beam> &outBeams) override;
	void SplitBeamByParticle(double, double, const std::vector<std::vector<int>> &,
							 std::vector<Beam> &) override; ///> for predefined trajectories

protected:
	void TraceInternalBeams(std::vector<Beam> &outBeams);
};
