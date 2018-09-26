#pragma once

#include "TracerPO.h"

class TracerBackScatterPoint : public TracerPO
{
public:
	TracerBackScatterPoint(Particle *particle, int reflNum,
						   const std::string &resultFileName);

	void Trace(const AngleRange &betaRange, const AngleRange &gammaRange,
			   const Tracks &tracks, double wave);

private:
	std::string GetTableHead(const AngleRange &range);
	void CreateResultFiles(ScatteringFiles &files, const Tracks &tracks,
						const std::string &subdir, const std::string &prefix = "");
	void CreateGroupResultFiles(const Tracks &tracks, ScatteringFiles &files,
								const std::string &subdir,
								const std::string &prefix = "");
	void AllocGroupMatrices(std::vector<Arr2D> &mtrcs, int maxGroupID);

	bool isNan = false;
};
