#include "HandlerTotalGO.h"

<<<<<<< HEAD

HandlerTotalGO::HandlerTotalGO(Particle *particle, Light *incidentLight,
							   float wavelength)
=======
using namespace std;

HandlerTotalGO::HandlerTotalGO(Particle *particle, Light *incidentLight, float wavelength)
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46
	: HandlerGO(particle, incidentLight, wavelength)
{
}

void HandlerTotalGO::HandleBeams(std::vector<Beam> &beams)
{
<<<<<<< HEAD
	m_sinZenith = sin(m_particle->rotAngle.beta);
=======
	m_sinAngle = sin(m_particle->rotAngle.zenith);
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46

	for (Beam &beam : beams)
	{
		beam.RotateSpherical(-m_incidentLight->direction,
							 m_incidentLight->polarizationBasis);
		// absorbtion
<<<<<<< HEAD
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
=======
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
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46
