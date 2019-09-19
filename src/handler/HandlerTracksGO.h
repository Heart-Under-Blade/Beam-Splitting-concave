#pragma once

#include "HandlerGO.h"

class HandlerTracksGO : public HandlerGO
{
public:
	HandlerTracksGO(Scattering *scattering, double wavelength = 0);

	void HandleBeams(std::vector<Beam> &beams) override;
	void WriteMatricesToFile(std::string &destName) override;
};
