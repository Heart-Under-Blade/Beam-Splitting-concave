#include "HandlerTracksGO.h"

HandlerTracksGO::HandlerTracksGO(Scattering *scattering,
								 double wavelength)
	: HandlerGO(scattering, wavelength)
{
}

void HandlerTracksGO::HandleBeams(std::vector<Beam> &beams)
{
	m_sinZenith = sin(m_particle->rotAngle.zenith);

	for (Beam &beam : beams)
	{
		int groupId = m_tracks->FindGroupByTrackId(beam.id);

		if (groupId >= 0)
		{
			beam.RotateSpherical(-m_startBeam.direction,
								 m_startBeam.polarizationBasis);

			const float &z = beam.direction.coordinates[2];
			int zenith = round(Orientation::RadToDeg(acos(z)));
			matrix m = ComputeMueller(zenith, beam);

			m_totalContrib.AddMueller(z, zenith, m);
			m_tracksContrib[groupId].AddMueller(z, zenith, m);
		}
	}
}

void HandlerTracksGO::WriteMatricesToFile(std::string &destName)
{
//	string dir = CreateFolder(destName);
//	dir += destName + "\\";

	for (size_t i = 0; i < m_tracksContrib.size(); ++i)
	{
		if ((*m_tracks)[i].size != 0)
		{
			std::string subname = (*m_tracks)[i].CreateGroupName();
//			AverageOverAlpha(true, m_normIndex, m_tracksContrib[i]);
			WriteToFile(m_tracksContrib[i], m_normIndex, destName + '_' +  subname);
		}
	}

	AverageOverAlpha(true, m_normIndex, m_totalContrib, destName);
	WriteToFile(m_totalContrib, m_normIndex, destName + "_all");
}
