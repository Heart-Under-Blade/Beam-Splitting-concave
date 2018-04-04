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
	void OutputContribution(const Tracks &tracks, PointContribution &contrib,
							ScatteringFiles &files, double degree, std::string prefix = "");
	void CreateResultFiles(ScatteringFiles &files, const Tracks &tracks,
						const std::string &subdir, const std::string &prefix = "");
	void CreateGroupResultFiles(const Tracks &tracks, ScatteringFiles &files,
								const std::string &subdir,
								const std::string &prefix = "");
	void AllocGroupMatrices(std::vector<Arr2D> &mtrcs, size_t maxGroupID);

	bool isNan = false;
};
