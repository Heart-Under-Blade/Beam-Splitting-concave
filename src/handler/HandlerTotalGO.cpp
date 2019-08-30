#include "HandlerTotalGO.h"

HandlerTotalGO::HandlerTotalGO(Particle *particle, Light *incidentLight,
							   float wavelength)
	: HandlerGO(particle, incidentLight, wavelength)
{
}

void HandlerTotalGO::HandleBeams(std::vector<Beam> &beams)
{
	m_sinZenith = sin(m_particle->rotAngle.zenith);

	for (Beam &beam : beams)
	{
		beam.RotateSpherical(-m_incidentLight->direction,
							 m_incidentLight->polarizationBasis);
		// absorption
		if (m_hasAbsorption && beam.actNo > 0)
		{
			ApplyAbsorption(beam);
		}

		const float &z = beam.direction.coordinates[2];
		int zenith = round(Orientation::RadToDeg(acos(z)));
		matrix m = ComputeMueller(zenith, beam);
		m_totalContrib.AddMueller(z, zenith, m);
	}
}

void HandlerTotalGO::WriteMatricesToFile(std::string &destName)
{
	AverageOverAlpha(true, m_normIndex, m_totalContrib, destName);
	WriteToFile(m_totalContrib, m_normIndex, destName + "_all");
}
