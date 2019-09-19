#pragma once

#include "HandlerGO.h"

class HandlerTotalGO : public HandlerGO
{
public:
	HandlerTotalGO(Scattering *scattering, double wavelength = 0);

	void HandleBeams(std::vector<Beam> &beams) override;
	void WriteMatricesToFile(std::string &destName) override;
};
