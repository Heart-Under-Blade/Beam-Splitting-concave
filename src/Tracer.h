#pragma once

#include "geometry_lib.h"
#include "CalcTimer.h"
#include "Mueller.hpp"
#include "BigInteger.hh"
#include "Handler.h"

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

class Tracer
{
public:
	Tracer(Particle *particle, int nActs, const std::string &resultFileName);
	~Tracer();

	// REF: delete?
	void TraceRandomPO2(int betaNumber, int gammaNumber, const ScatteringSphere &bsCone,
						const Tracks &tracks, double wave);

	void SetHandler(Handler *handler);

	void SetIsOutputGroups(bool value);// REF: заменить

	void OutputStatisticsPO(CalcTimer &timer, long long orNumber, const std::string &path);

	Light m_incidentLight;
	std::string m_summary;

protected:
	Handler *m_handler;
	Scattering *m_scattering;
	Particle *m_particle;

	double m_incomingEnergy;
	double m_outcomingEnergy;

	std::string m_resultDirName;
	double m_wavelength;
	Symmetry m_symmetry;
	time_t m_startTime;

	// REF: заменить
	bool isOutputGroups = false;

	long long m_timeElapsed = 0;

protected:
	void OutputStartTime(CalcTimer &timer);
	void OutputProgress(int nOrientation, long long count,
						int zenith, int azimuth, CalcTimer &timer);
	void OutputOrientationToLog(int i, int j, std::ostream &logfile);

private:
	void HandleBeamsPO2(std::vector<Beam> &outBeams, const ScatteringSphere &bsCone, int groupID);
	void SetIncidentLight(Particle *particle);
};
