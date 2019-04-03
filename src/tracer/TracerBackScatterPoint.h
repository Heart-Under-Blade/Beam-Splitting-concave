#pragma once

#include "TracerPO.h"

class TracerBackScatterPoint : public TracerPO
{
public:
	TracerBackScatterPoint(Particle *particle, Scattering *scattering,
						   const std::string &resultFileName);

	void TraceRandom(const AngleRange &betaRange, const AngleRange &gammaRange) override;

private:
	std::string GetTableHead(const AngleRange &range);
	void CreateResultFiles(ScatteringFiles &files, Tracks *tracks,
						   const std::string &subdir, const std::string &prefix = "");
	void CreateGroupResultFiles(ScatteringFiles &files, Tracks *tracks,
								const std::string &subdir,
								const std::string &prefix = "");
	void AllocGroupMatrices(std::vector<Arr2D> &mtrcs, int maxGroupID);

	bool isNan = false;
};
