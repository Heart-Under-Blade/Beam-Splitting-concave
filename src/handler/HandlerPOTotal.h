#pragma once

#include "HandlerPO.h"

class HandlerPOTotal : public HandlerPO
{
public:
	HandlerPOTotal(Particle *particle, Light *incidentLight,
				   double wavelength = 0);
    void WriteMatricesToFile(std::string &destName) override;
	void AddToMueller() override;

    matrix *m_Lp;
    matrix *m_Ln;
};
