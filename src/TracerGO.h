#pragma once

#include "Tracer.h"
#include "Handler.h"

class TracerGO : public Tracer
{
public:
	TracerGO(Particle *particle, int reflNum, const std::string &resultFileName);

	void TraceRandom(const AngleRange &betaRange, const AngleRange &gammaRange,
					 bool isCalcTracks, bool isAbs, double wavelenght = 0);

	void TraceFixed(const double &beta, const double &gamma,
					bool isCalcTracks, bool isAbs, double wavelenght = 0);
	void SetTracks(Tracks *tracks);

protected:
	Tracks *m_tracks;
	HandlerGO *m_handler;

protected:
	double CalcNorm(long long orNum);
	void OutputSummary(int orNumber, double D_tot, double NRM, CalcTimer &timer);

private:
	void SetHandler(bool isCalcTracks, bool isAbs, double wavelength = 0);
};
