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
	Tracer(Particle *particle, int maxActNo, const std::string &resultFileName);
	~Tracer();

	void SetHandler(Handler *handler);

	void SetIsOutputGroups(bool value);// REF: заменить

	void OutputStatisticsPO(CalcTimer &timer, long long orNumber, const std::string &path);

	Light m_incidentLight;

protected:
	Particle *m_particle;
	Handler *m_handler;
	Scattering *m_scattering;

	double m_incomingEnergy;
	double m_outcomingEnergy;

	std::string m_resultDirName;
	double m_wavelength;
	std::string m_summary;
	time_t m_startTime;

	// REF: заменить
	bool isOutputGroups = false;

protected:
	void OutputStartTime(CalcTimer &timer);
	void OutputProgress(int betaNumber, long long count, CalcTimer &timer);
	void OutputOrientationToLog(int i, int j, std::ostream &logfile);

private:
	void SetIncidentLight(Particle *particle);
};
