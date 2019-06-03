#include "HandlerTotalGO.h"


HandlerTotalGO::HandlerTotalGO(Particle *particle, Light *incidentLight,
							   float wavelength)
	: HandlerGO(particle, incidentLight, wavelength)
{
}

void HandlerTotalGO::HandleBeams(std::vector<Beam> &beams)
{
	m_sinZenith = sin(m_particle->rotAngle.beta);

	for (Beam &beam : beams)
	{
		beam.RotateSpherical(-m_incidentLight->direction,
							 m_incidentLight->polarizationBasis);
		// absorbtion
		if (m_hasAbsorption && beam.nActs > 0)
		{
			ApplyAbsorption(beam);
		}

		const float &z = (acos(beam.direction.cz)*SPHERE_RING_NUM)/M_PI;
		int zenith = round(z);
		matrix m = ComputeMueller(zenith, beam);

		m_totalContrib.AddMueller(z, zenith, m);
	}
}

void HandlerTotalGO::WriteMatricesToFile(std::string &destName)
{
	AverageOverAlpha(true, m_normIndex, m_totalContrib, destName);
	WriteToFile(m_totalContrib, m_normIndex, destName + "_all");
}
