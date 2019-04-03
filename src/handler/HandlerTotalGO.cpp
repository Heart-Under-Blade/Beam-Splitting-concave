#include "HandlerTotalGO.h"

using namespace std;

HandlerTotalGO::HandlerTotalGO(Particle *particle, Light *incidentLight, float wavelength)
	: HandlerGO(particle, incidentLight, wavelength)
{
}

void HandlerTotalGO::HandleBeams(std::vector<Beam> &beams)
{
	m_sinAngle = sin(m_particle->rotAngle.zenith);

	for (Beam &beam : beams)
	{
		beam.RotateSpherical(-m_incidentLight->direction,
							 m_incidentLight->polarizationBasis);
		// absorbtion
		if (m_hasAbsorbtion && beam.actNo > 0)
		{
			ApplyAbsorbtion(beam);
		}

		const float &z = beam.direction.coordinates[2];
		int zenith = round(Orientation::RadToDeg(acos(z)));
		matrix m = ComputeMueller(zenith, beam);
		m_totalContrib.AddMueller(z, zenith, m);
	}
}
