#include "HandlerTracksGO.h"


HandlerTracksGO::HandlerTracksGO(Particle *particle, Light *incidentLight, float wavelength)
	: HandlerGO(particle, incidentLight, wavelength)
{
}

void HandlerTracksGO::HandleBeams(std::vector<Beam> &beams)
{
	m_sinZenith = sin(m_particle->rotAngle.beta);

	for (Beam &beam : beams)
	{
#ifdef _DEBUG // DEB
//        std::vector<int> tr;
//        Tracks::RecoverTrack(beam, m_particle->nFacets, tr);
#endif
		int groupId = m_tracks->FindGroupByTrackId(beam.id);

		if (groupId >= 0 || !m_tracks->shouldComputeTracksOnly)
		{
			beam.RotateSpherical(-m_incidentLight->direction,
								 m_incidentLight->polarizationBasis);
			// absorbtion
			if (m_hasAbsorption && beam.nActs > 0)
			{
				ApplyAbsorption(beam);
			}

			const float &z = beam.direction.cz;
			int zenith = round((acos(z)*SPHERE_RING_NUM)/M_PI);
			matrix m = ComputeMueller(z, beam);

			m_totalContrib->AddMueller(z, zenith, m);
			m_tracksContrib[groupId]->AddMueller(z, zenith, m);
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
			WriteToFile(*m_tracksContrib[i], m_normIndex, destName + '_' +  subname);
		}
	}

	AverageOverAlpha(true, m_normIndex, *m_totalContrib, destName);
	WriteToFile(*m_totalContrib, m_normIndex, destName + "_all");
}



