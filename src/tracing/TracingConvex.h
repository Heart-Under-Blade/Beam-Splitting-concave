#pragma once

#include "Tracing.h"

class TracingConvex : public Tracing
{
public:
	TracingConvex(Particle *particle, const Point3f &incidentBeamDir,
				  bool isOpticalPath, const Point3f &polarizationBasis,
				  int interReflectionNumber);

	void SplitBeamByParticle(std::vector<OutBeam> &outBeams,
							 double &lightSurfaceSquare) override;

	void SplitBeamByParticle(const std::vector<std::vector<int>> &tracks,
							 std::vector<OutBeam> &outBeams) override; ///> for predefined trajectories

};
