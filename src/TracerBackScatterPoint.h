#pragma once

#include "Tracer.h"

struct AngleRange
{
	double min;
	double max;
	int number;
	double norm;
	double step;

	AngleRange(double _min, double _max, int _number)
		: number(_number)
	{
		min = _min;
		max = _max;
		norm = max - min;
		step = norm/number;
	}
};

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
								const std::string &prefix = "");
	void AllocGroupMatrices(std::vector<Arr2D> &mtrcs, size_t maxGroupID);
};
