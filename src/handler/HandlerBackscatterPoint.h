#pragma once

#include "HandlerPO.h"

class HandlerBackScatterPoint : public HandlerPO
{
public:
	HandlerBackScatterPoint(Particle *particle, Light *incidentLight,
							double wavelength = 0);

	void HandleBeams(std::vector<Beam> &beams) override;
	void SetTracks(Tracks *tracks) override;

	void OutputContribution(ScatteringFiles &files, double angle, double energy,
							bool isOutputGroups, std::string prefix = "");
private:
	PointContribution *originContrib;
	PointContribution *correctedContrib;
};

