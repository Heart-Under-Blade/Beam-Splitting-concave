#include "HandlerTotalGO.h"

HandlerTotalGO::HandlerTotalGO(Scattering *scattering, double wavelength)
	: HandlerGO(scattering, wavelength)
{
}

void HandlerTotalGO::HandleBeams(std::vector<Beam> &beams)
{
	m_sinZenith = sin(m_particle->rotAngle.zenith);

	for (Beam &beam : beams)
	{
		beam.RotateSpherical(-m_startBeam.direction,
							 m_startBeam.polarizationBasis);

		// absorption
		if (m_hasAbsorption && beam.actNo > 0)
		{
			BeamInfo info = ComputeBeamInfo(beam);
			ApplyAbsorption(info);
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
