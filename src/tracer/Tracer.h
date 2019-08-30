#pragma once

#include "geometry_lib.h"
#include "CalcTimer.h"
#include "Mueller.hpp"
#include "BigInteger.hh"
#include "ArgumentParser.h"
#include "Scattering.h"
#include "Handler.h"

/**
 * @brief Scatters the light on a particle, collect result beams and handle them
 */
class Tracer
{
public:
	Tracer(Particle *particle, Scattering *scattering,
				const std::string &resultFileName);
	~Tracer();

	/**
	 * @brief Trace a light on a random oriented particle rotated by given angle ranges
	 */
	virtual void TraceRandom(const OrientationRange &/*range*/);
	/**
	 * @brief Trace a light on a fixed orienteted particle with given orientation
	 * @param orientation value of a particle orientation
	 */
	void TraceFixed(const Orientation &orientation);

	void SetHandler(Handler *handler);

	void SetIsOutputGroups(bool value);// REF: заменить

	void OutputStatisticsPO(CalcTimer &timer, long long orNumber, const std::string &path);

public:
	Light m_incidentLight;
	std::string m_summary;

protected:
	Particle *m_particle;
	Handler *m_handler;
	Scattering *m_scattering;

	double m_incomingEnergy;
	double m_outcomingEnergy;

	std::string m_resultDirName;
	double m_wavelength;
	Orientation m_symmetry;
	time_t m_startTime;

	// REF: заменить
	bool isOutputGroups = false;

	long long m_timeElapsed = 0;

protected:
	void OutputStartTime(CalcTimer &timer);
	void OutputProgress(long long nOrientation, long long count,
						CalcTimer &timer);
	void OutputOrientationToLog(int i, int j, std::ostream &logfile);

private:
};
