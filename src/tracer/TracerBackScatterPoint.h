#pragma once

#include "TracerPO.h"

class TracerBackScatterPoint : public TracerPO
{
public:
	TracerBackScatterPoint(Particle *particle, Scattering *scattering,
						   const std::string &resultFileName);

	void TraceRandom(const OrientationRange &range) override;

private:
	std::string GetTableHead(const OrientationRange &range);
	void CreateResultFiles(ScatteringFiles &files, Tracks *tracks,
						   const std::string &subdir,
						   const std::string &prefix = "");
	void CreateGroupResultFiles(ScatteringFiles &files, Tracks *tracks,
								const std::string &subdir,
								const std::string &prefix = "");
	void AllocGroupMatrices(std::vector<Arr2D> &mtrcs, size_t maxGroupID);

	bool isNan = false;
};
