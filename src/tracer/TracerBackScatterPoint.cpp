#include "TracerBackScatterPoint.h"

#include "common.h"
#include "ScatteringFiles.h"
#include "HandlerBackscatterPoint.h"

#include <iostream>

using namespace std;

TracerBackScatterPoint::TracerBackScatterPoint(Particle *particle, Scattering *scattering,
											   const std::string &resultFileName)
	: TracerPO(particle, scattering, resultFileName)
{
}

void TracerBackScatterPoint::TraceRandom(const AngleRange &betaRange,
										 const AngleRange &gammaRange)
{
	int nGroups = m_handler->m_tracks->size();

	m_wavelength = m_handler->m_wavelength;
	CalcTimer timer;
	long long count = 0;

	string dir = CreateFolder(m_resultDirName);
	string tableHead = GetTableHead(betaRange);

	string fulldir = dir + m_resultDirName + '\\';
	ofstream logfile(fulldir + "log.txt", ios::out);

	ScatteringFiles resFiles(dir, m_resultDirName, tableHead);
	CreateDir(fulldir + "res");
	CreateResultFiles(resFiles, m_handler->m_tracks, "res");

	vector<Arr2D> groupResultM;
	AllocGroupMatrices(groupResultM, nGroups);

	ScatteringFiles corFiles(dir, m_resultDirName, tableHead);
	CreateDir(fulldir + "cor");
	CreateResultFiles(corFiles, m_handler->m_tracks, "cor", "cor_");

	vector<Arr2D> groupResultM_cor;
	AllocGroupMatrices(groupResultM_cor, nGroups);

	vector<Beam> outBeams;

	OutputStartTime(timer);

	Orientation angle;

	for (int i = 0; i <= betaRange.number; ++i)
	{
		m_incomingEnergy = 0;
		OutputProgress(betaRange.number, ++count, timer);

		angle.zenith = betaRange.min + betaRange.step*i;

		for (int j = 0; j <= gammaRange.number; ++j)
		{
			angle.azimuth = gammaRange.min + gammaRange.step*j;
#ifdef _DEBUG // DEB
//			angle.beta = Angle::DegToRad(179.34);
//			angle.gamma = Angle::DegToRad(37);
#endif
			m_particle->Rotate(angle);
			m_scattering->ScatterLight(outBeams);
#ifdef _DEBUG // DEB
			m_particle->Output();
#endif
			m_incomingEnergy += m_scattering->GetIncidentEnergy();

			m_handler->HandleBeams(outBeams);
			outBeams.clear();
			OutputOrientationToLog(i, j, logfile);		

			if (isNan)
			{
				logfile << "------------------";
				isNan = false;
				continue;
			}
		}

		double degBeta = Orientation::RadToDeg(angle.zenith);

		// REF: remove static casts
		static_cast<HandlerBackScatterPoint*>
				(m_handler)->OutputContribution(resFiles, degBeta,
												m_incomingEnergy, isOutputGroups);
		static_cast<HandlerBackScatterPoint*>
				(m_handler)->OutputContribution(corFiles, degBeta,
												m_incomingEnergy, isOutputGroups,
												"cor_");
	}

	EraseConsoleLine(50);
	cout << "100%";

	if (isOutputGroups)
	{
		for (int group = 0; group < nGroups; ++group)
		{
			ofstream &file = *(resFiles.GetGroupFile(group));
			file.close();

			ofstream &file_cor = *(corFiles.GetGroupFile(group));
			file_cor.close();
		}
	}

	long long orNumber = betaRange.number * gammaRange.number;
	OutputStatisticsPO(timer, orNumber, m_resultDirName);
}

string TracerBackScatterPoint::GetTableHead(const AngleRange &range)
{
	return to_string(range.number) + ' '
			+ to_string(Orientation::RadToDeg(range.max)) + ' '
			+ to_string(Orientation::RadToDeg(range.step)) + '\n'
			+ "beta cr_sec M11 M12 M13 M14 M21 M22 M23 M24 M31 M32 M33 M34 M41 M42 M43 M44"
			+ '\n';
}

void TracerBackScatterPoint::CreateResultFiles(ScatteringFiles &files,
											   Tracks *tracks,
											   const string &subdir,
											   const string &prefix)
{
	files.CreateMainFile(subdir, prefix + "all");
	files.CreateMainFile(subdir, prefix + "other");
	files.CreateMainFile(subdir, prefix + "difference");

	if (isOutputGroups)
	{
		CreateGroupResultFiles(files, tracks, subdir, prefix);
	}
}

void TracerBackScatterPoint::CreateGroupResultFiles(ScatteringFiles &files,
													Tracks *tracks,
													const string &subdir,
													const string &prefix)
{
	for (int i = 0; i < tracks->size(); ++i)
	{
		string groupName = (*tracks)[i].CreateGroupName();
		string filename = prefix + groupName;
		files.CreateGroupFile(subdir, filename);
	}
}

void TracerBackScatterPoint::AllocGroupMatrices(vector<Arr2D> &mtrcs, int maxGroupID)
{
	for (int i = 0; i < maxGroupID; ++i)
	{
		Arr2D m(1, 1, 4, 4);
		mtrcs.push_back(m);
	}
}
