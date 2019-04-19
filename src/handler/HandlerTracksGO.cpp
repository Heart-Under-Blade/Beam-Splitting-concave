#include "HandlerTracksGO.h"

using namespace std;

HandlerTracksGO::HandlerTracksGO(Particle *particle, Light *incidentLight, float wavelength)
	: HandlerGO(particle, incidentLight, wavelength)
{
}

void HandlerTracksGO::HandleBeams(std::vector<Beam> &beams)
{
	m_sinAngle = sin(m_particle->rotAngle.zenith);

	for (Beam &beam : beams)
	{
		int groupId = m_tracks->FindGroupByTrackId(beam.id);

		if (groupId >= 0)
		{
			beam.RotateSpherical(-m_incidentLight->direction,
								 m_incidentLight->polarizationBasis);

			const float &z = beam.direction.coordinates[2];
			int zenith = round(Orientation::RadToDeg(acos(z)));
			matrix m = ComputeMueller(zenith, beam);

			m_totalContrib.AddMueller(z, zenith, m);
			m_tracksContrib[groupId].AddMueller(z, zenith, m);
		}
	}
}

void HandlerTracksGO::WriteMatricesToFile(string &destName)
{
	for (size_t i = 0; i < m_tracksContrib.size(); ++i)
	{
		if ((*m_tracks)[i].size != 0)
		{
			string subname = (*m_tracks)[i].CreateGroupName();
			WriteToFile(m_tracksContrib[i], m_normIndex, destName + '_' +  subname);
		}
	}

	HandlerGO::WriteMatricesToFile(destName);
}
