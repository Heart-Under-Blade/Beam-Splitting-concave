#pragma once

class Beam;
class Splitting;

class Incidence
{
public:
	virtual void ComputeLightParams(const Beam &incidentBeam,
									Splitting &splitter) const {};

	virtual void ComputeJonesMatrices(const Beam &incidentBeam,
									Splitting &splitter) const {};

	void ComputeOpticalPaths(const Beam &incidentBeam, Splitting &splitter) const;
};
