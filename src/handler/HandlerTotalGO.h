#pragma once

#include "HandlerGO.h"

class HandlerTotalGO : public HandlerGO
{
public:
	HandlerTotalGO(Particle *particle, Light *incidentLight, float wavelength = 0);

	void HandleBeams(std::vector<Beam> &beams) override;
<<<<<<< HEAD
	void WriteMatricesToFile(std::string &destName) override;
};

=======
};
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46
