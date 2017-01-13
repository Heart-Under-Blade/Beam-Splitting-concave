#pragma once

#include "Tracing.h"

class TracingConvex : public Tracing
{
public:
	TracingConvex(Particle *particle, const Point3f &incidentBeamDir,
				  bool isOpticalPath, const Point3f &polarizationBasis,
				  int interReflectionNumber);

	void SplitBeamByParticle(std::vector<Beam> &outBeams) override;

	double BeamCrossSection(const Beam &beam) const override;

	void SplitBeamByParticle(const std::vector<std::vector<int>> &,
							 std::vector<Beam> &) override; ///> for predefined trajectories

protected:
	void TraceInternalReflections(std::vector<Beam> &outBeams);
};
