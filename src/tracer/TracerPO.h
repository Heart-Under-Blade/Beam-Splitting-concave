#pragma once

#include "Tracer.h"

class TracerPO : public LightTracer
{
public:
<<<<<<< HEAD
	TracerPO(Particle *particle, int nActs, const std::string &resultFileName);
	virtual void TraceRandom(const AngleRange &betaRange, const AngleRange &gammaRange);
	void TraceFixed(const double &beta, const double &gamma);
=======
	TracerPO(Particle *particle, Scattering *scattering,
			 const std::string &resultFileName);

	void TraceRandom(const AngleRange &zenithRange,
					 const AngleRange &azimuthRange);
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46
};
