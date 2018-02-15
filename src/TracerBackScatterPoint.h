#pragma once

#include "Tracer.h"

class TracerBackScatterPoint : public Tracer
{
public:
	TracerBackScatterPoint(Particle *particle, int reflNum,
						   const std::string &resultFileName);

	void Trace(const AngleRange &betaRange, const AngleRange &gammaRange,
			   const Tracks &tracks, double wave);

protected:
	void HandleBeams(std::vector<Beam> &outBeams, const Tracks &tracks,
					 PointContribution &general, PointContribution &corrected);

private:
	std::string GetTableHead(const AngleRange &range);
	void OutputContribution(size_t groupNumber, PointContribution &contrib,
							ScatteringFiles &files, double degree, std::string prefix = "");
	void CreateResultFiles(ScatteringFiles &files, const Tracks &tracks,
						const std::string &subdir, const std::string &prefix = "");
	void CreateGroupResultFiles(const Tracks &tracks, ScatteringFiles &files,
								const std::string &subdir,
								const std::string &prefix = "");
	void AllocGroupMatrices(std::vector<Arr2D> &mtrcs, size_t maxGroupID);
};
