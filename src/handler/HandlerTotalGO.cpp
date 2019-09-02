#include "HandlerTotalGO.h"

<<<<<<< HEAD

=======
>>>>>>> origin/refactor
HandlerTotalGO::HandlerTotalGO(Particle *particle, Light *incidentLight,
							   float wavelength)
	: HandlerGO(particle, incidentLight, wavelength)
{
}

void HandlerTotalGO::HandleBeams(std::vector<Beam> &beams)
{
<<<<<<< HEAD
	m_sinZenith = sin(m_particle->rotAngle.beta);
=======
	m_sinZenith = sin(m_particle->rotAngle.zenith);
>>>>>>> origin/refactor

	for (Beam &beam : beams)
	{
		beam.RotateSpherical(-m_incidentLight->direction,
							 m_incidentLight->polarizationBasis);
<<<<<<< HEAD
		// absorbtion
		if (m_hasAbsorption && beam.nActs > 0)
=======
		// absorption
		if (m_hasAbsorption && beam.actNo > 0)
>>>>>>> origin/refactor
		{
			ApplyAbsorption(beam);
		}

<<<<<<< HEAD
		const float &z = (acos(beam.direction.cz)*SPHERE_RING_NUM)/M_PI;
		int zenith = round(z);
		matrix m = ComputeMueller(zenith, beam);

=======
		const float &z = beam.direction.coordinates[2];
		int zenith = round(Orientation::RadToDeg(acos(z)));
		matrix m = ComputeMueller(zenith, beam);
>>>>>>> origin/refactor
		m_totalContrib.AddMueller(z, zenith, m);
	}
}

void HandlerTotalGO::WriteMatricesToFile(std::string &destName)
{
	AverageOverAlpha(true, m_normIndex, m_totalContrib, destName);
	WriteToFile(m_totalContrib, m_normIndex, destName + "_all");
}
