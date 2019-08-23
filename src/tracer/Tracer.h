#pragma once

#include "geometry_lib.h"
#include "CalcTimer.h"
#include "Mueller.hpp"
#include "BigInteger.hh"
#include "Handler.h"
#include "ArgumentParser.h"

/**
 * @brief Scatters the light on a particle, collect result beams and handle them
 */
class LightTracer
{
public:
<<<<<<< HEAD
	Tracer(Particle *particle, int nActs, const std::string &resultFileName);
	~Tracer();

	// REF: delete?
//	void TraceRandomPO2(int betaNumber, int gammaNumber, const Conus &bsCone,
//						const Tracks &tracks, double wave);
=======
	LightTracer(Particle *particle, Scattering *scattering,
				const std::string &resultFileName);
	~LightTracer();

	/**
	 * @brief Trace a light on a random oriented particle rotated by given angle ranges
	 * @param zenithRange range for rotation of particle by zinith angle
	 * @param azimuthRange range for rotation of particle by azimuth angle
	 */
	virtual void TraceRandom(const AngleRange &zenithRange,
							 const AngleRange &azimuthRange);
	/**
	 * @brief Trace a light on a fixed orienteted particle with given orientation
	 * @param orientation value of a particle orientation
	 */
	void TraceFixed(const Orientation &orientation);
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46

	void SetHandler(Handler *handler);

	void SetIsOutputGroups(bool value);// REF: заменить

	void OutputLogPO(CalcTimer &timer, long long orNumber, const std::string &path);

<<<<<<< HEAD
	Light m_incidentLight;
	std::string m_log;

=======
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46
protected:
	Particle *m_particle;
	Handler *m_handler;
	Scattering *m_scattering;

	double m_incomingEnergy;
	double m_outcomingEnergy;

	std::string m_resultDirName;
	double m_wavelength;
<<<<<<< HEAD
	Symmetry m_symmetry;
=======
	std::string m_summary;
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46
	time_t m_startTime;

	// REF: заменить
	bool isOutputGroups = false;

protected:
	void OutputStartTime(CalcTimer &timer);
	void OutputProgress(int betaNumber, long long count, CalcTimer &timer);
	void OutputOrientationToLog(int i, int j, std::ostream &logfile);

private:
<<<<<<< HEAD
//	void HandleBeamsPO2(std::vector<Beam> &outBeams, const Conus &bsCone, int groupID);
=======
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46
	void SetIncidentLight(Particle *particle);
};
