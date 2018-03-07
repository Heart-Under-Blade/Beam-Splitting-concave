#pragma once

#include "Tracer.h"
#include "Handler.h"

class TracerGO : public Tracer
{
public:
	TracerGO(Particle *particle, int reflNum, const std::string &resultFileName);

	void TraceRandom(const AngleRange &betaRange, const AngleRange &gammaRange);
	void TraceFixed(const double &beta, const double &gamma);

	void SetHandler(HandlerGO *handler);

protected:
	HandlerGO *m_handler;

protected:
	double CalcNorm(long long orNum);
	void OutputSummary(int orNumber, double D_tot, double NRM, CalcTimer &timer);
};
