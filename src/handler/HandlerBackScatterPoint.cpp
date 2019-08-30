#include "HandlerBackScatterPoint.h"

#define BEAM_DIR_LIM 0.9396 // cos(20)

HandlerBackScatterPoint::HandlerBackScatterPoint(Particle *particle,
												 Light *incidentLight,
												 double wavelength)
	: HandlerPO(particle, incidentLight, wavelength)
{
}

void HandlerBackScatterPoint::HandleBeams(std::vector<Beam> &beams)
{
	Point3d backDirection(0, 0, 1);
	Point3d vf = -m_incidentLight->polarizationBasis;

	for (Beam &beam : beams)
	{
		if (beam.direction.cz < BEAM_DIR_LIM)
		{
			continue;
		}

		int groupId = m_tracks->FindGroupByTrackId(beam.id);

		if (groupId < 0 && m_tracks->shouldComputeTracksOnly)
		{
			continue;
		}

		beam.polarizationBasis = beam.RotateSpherical(-m_incidentLight->direction,
													  m_incidentLight->polarizationBasis);

		BeamInfo info = ComputeBeamInfo(beam);

		matrixC diffractedMatrix = ApplyDiffraction(beam, info, backDirection, vf);

		// correction
		Matrix2x2c jonesCor = diffractedMatrix;
		jonesCor.m12 -= jonesCor.m21;
		jonesCor.m12 /= 2;
		jonesCor.m21 = -jonesCor.m12;

		if (groupId < 0 && !m_tracks->shouldComputeTracksOnly)
		{
			originContrib->AddToMueller(diffractedMatrix);
			correctedContrib->AddToMueller(jonesCor);
		}
		else
		{
			originContrib->AddToGroup(diffractedMatrix, groupId);
			correctedContrib->AddToGroup(jonesCor, groupId);
		}
	}

	originContrib->SumGroupTotal();
	correctedContrib->SumGroupTotal();
}

void HandlerBackScatterPoint::SetTracks(Tracks *tracks)
{
	Handler::SetTracks(tracks);
	originContrib = new PointContribution(tracks->size(), m_normIndex);
	correctedContrib = new PointContribution(tracks->size(), m_normIndex);
}

void HandlerBackScatterPoint::OutputContribution(ScatteringFiles &files,
												 double angle, double energy,
												 bool isOutputGroups,
												 std::string prefix)
{
	PointContribution *contrib = (prefix == "") ? originContrib : correctedContrib;
	contrib->SumTotal();

	energy *= m_normIndex;

	std::ofstream *all = files.GetMainFile(prefix + "all");
	*(all) << angle << ' ' << energy << ' ';
	*(all) << contrib->GetTotal() << std::endl;
//cout << endl << endl << contrib->GetRest()(0,0) << endl << endl ;
	if (isOutputGroups)
	{
		for (size_t gr = 0; gr < m_tracks->size(); ++gr)
		{
			std::ofstream &file = *(files.GetGroupFile(gr));
			file << angle << ' ' << energy << ' ';
			file << contrib->GetGroupMueller(gr) << std::endl;
		}
	}

	if (!m_tracks->shouldComputeTracksOnly)
	{
		std::ofstream &other = *(files.GetMainFile(prefix + "other"));
		other << angle << ' ' << energy << ' ';
		other << contrib->GetRest() << std::endl;

		std::ofstream &diff = *(files.GetMainFile(prefix + "difference"));
		diff << angle << ' ' << energy << ' ';
		diff << contrib->GetGroupTotal() << std::endl;
	}

	contrib->Reset();
}
