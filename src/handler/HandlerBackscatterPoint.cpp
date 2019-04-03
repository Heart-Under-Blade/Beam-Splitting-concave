#include "HandlerBackscatterPoint.h"

#define BEAM_DIR_LIM		0.9396

using namespace std;

HandlerBackScatterPoint::HandlerBackScatterPoint(Particle *particle,
												 Light *incidentLight,
												 float wavelength)
	: HandlerPO(particle, incidentLight, wavelength)
{
}

void HandlerBackScatterPoint::HandleBeams(std::vector<Beam> &beams)
{
	Point3d vr(0, 0, 1);
	Point3d vf = -m_incidentLight->polarizationBasis;

#ifdef _DEBUG // DEB
	int c = 0;
#endif
	for (Beam &beam : beams)
	{
#ifdef _DEBUG // DEB
		++c;
//		vector<int> tr;
//		m_tracks->RecoverTrack(beam, tr);
#endif

		if (con20 && beam.direction.coordinates[2] < BEAM_DIR_LIM)
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
#ifdef _DEBUG // DEB
//		double cross = BeamCrossSection(beam);
//		double area = cross*m_sinAngle;
//		vector<int> track;
//		double ddd = m[0][0];
		int zenith = round((acos(beam.direction.coordinates[2])*SPHERE_RING_NUM)/M_PI);
#endif
		// absorbtion
		if (m_hasAbsorbtion && beam.actNo > 0)
		{
			ApplyAbsorbtion(beam);
		}

		Point3f beamBasis = Point3f::CrossProduct(beam.polarizationBasis, beam.direction);
		beamBasis = beamBasis/Point3f::Length(beamBasis);

		Point3f center = beam.Center();
		auto path = ComputeOpticalPath(beam);
		double projLenght = path.GetTotal() + Point3f::DotProduct(center, beam.direction);

		matrixC jones(2, 2);
		matrixC fnJones = ComputeFnJones(beam.Jones, center, vr, projLenght);
		ApplyDiffraction(beam, beamBasis, vf, vr, fnJones, jones);

		// correction
		Matrix2x2c jonesCor = jones;

#ifdef _DEBUG // DEB
//		correctedContrib->AddToMueller(jonesCor);
#endif
		jonesCor.m12 -= jonesCor.m21;
		jonesCor.m12 /= 2;
		jonesCor.m21 = -jonesCor.m12;

		if (groupId < 0 && !m_tracks->shouldComputeTracksOnly)
		{
			originContrib->AddToMueller(jones);
			correctedContrib->AddToMueller(jonesCor);
		}
		else
		{
			originContrib->AddToGroup(jones, groupId);
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
												 string prefix)
{
	PointContribution *contrib = (prefix == "") ? originContrib : correctedContrib;
	contrib->SumTotal();

	energy *= m_normIndex;

	ofstream *all = files.GetMainFile(prefix + "all");
	*(all) << angle << ' ' << energy << ' ';
	*(all) << contrib->GetTotal() << endl;
//cout << endl << endl << contrib->GetRest()(0,0) << endl << endl ;
	if (isOutputGroups)
	{
		for (int gr = 0; gr < m_tracks->size(); ++gr)
		{
			ofstream &file = *(files.GetGroupFile(gr));
			file << angle << ' ' << energy << ' ';
			file << contrib->GetGroupMueller(gr) << endl;
		}
	}

	if (!m_tracks->shouldComputeTracksOnly)
	{
		ofstream &other = *(files.GetMainFile(prefix + "other"));
		other << angle << ' ' << energy << ' ';
		other << contrib->GetRest() << endl;

		ofstream &diff = *(files.GetMainFile(prefix + "difference"));
		diff << angle << ' ' << energy << ' ';
		diff << contrib->GetGroupTotal() << endl;
	}

	contrib->Reset();
}
