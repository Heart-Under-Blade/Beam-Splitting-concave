#pragma once

#include "HandlerGO.h"

class HandlerTracksGO : public HandlerGO
{
public:
<<<<<<< HEAD
    HandlerTracksGO(Particle *particle, Light *incidentLight, float wavelength = 0);

    void HandleBeams(std::vector<Beam> &beams) override;
    void WriteMatricesToFile(std::string &destName) override;
=======
	HandlerTracksGO(Particle *particle, Light *incidentLight, float wavelength = 0);

	void HandleBeams(std::vector<Beam> &beams) override;
	void WriteMatricesToFile(std::string &destName) override;
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46
};
