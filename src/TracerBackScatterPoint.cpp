#include "TracerBackScatterPoint.h"
#include "global.h"
#include "ScatteringFiles.h"
#include <iostream>

using namespace std;

TracerBackScatterPoint::TracerBackScatterPoint(Particle *particle, int reflNum,
											   const std::string &resultFileName)
	: TracerPO(particle, reflNum, resultFileName)
{
}

void TracerBackScatterPoint::Trace(const AngleRange &betaRange, const AngleRange &gammaRange,
								   const Tracks &tracks, double wave)
{
	size_t nGroups = tracks.size();
	int normIndex = gammaRange.step/gammaRange.norm;
	m_handler->SetNormIndex(normIndex);

	m_wavelength = wave;
	CalcTimer timer;
	long long count = 0;

	string dir = CreateFolder(m_resultDirName);
	string tableHead = GetTableHead(betaRange);

	string fulldir = dir + m_resultDirName + '\\';
	ofstream logfile(fulldir + "log.txt", ios::out);

	ScatteringFiles resFiles(dir, m_resultDirName, tableHead);
	CreateDir(fulldir + "res");
	CreateResultFiles(resFiles, tracks, "res");

	vector<Arr2D> groupResultM;
	AllocGroupMatrices(groupResultM, nGroups);

	ScatteringFiles corFiles(dir, m_resultDirName, tableHead);
	CreateDir(fulldir + "cor");
	CreateResultFiles(corFiles, tracks, "cor", "cor_");

	vector<Arr2D> groupResultM_cor;
	AllocGroupMatrices(groupResultM_cor, nGroups);

	vector<Beam> outBeams;

	OutputStartTime(timer);

	double beta, gamma;

	for (int i = 0; i <= betaRange.number; ++i)
	{
		m_incomingEnergy = 0;
		OutputProgress(betaRange.number, ++count, timer);

		beta = betaRange.min + betaRange.step*i;

		for (int j = 0; j <= gammaRange.number; ++j)
		{
			gamma = gammaRange.min + gammaRange.step*j;
			m_scattering->ScatterLight(beta, gamma, outBeams);

			m_incomingEnergy += m_scattering->GetIncomingEnergy();

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

		m_incomingEnergy *= normIndex;

		double degBeta = RadToDeg(beta);

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
		for (size_t group = 0; group < nGroups; ++group)
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
			+ to_string(RadToDeg(range.max)) + ' '
			+ to_string(RadToDeg(range.step)) + '\n'
			+ "beta cr_sec M11 M12 M13 M14 M21 M22 M23 M24 M31 M32 M33 M34 M41 M42 M43 M44"
			+ '\n';
}

void TracerBackScatterPoint::CreateResultFiles(ScatteringFiles &files, const Tracks &tracks,
											   const string &subdir,
											   const string &prefix)
{
	files.CreateMainFile(subdir, prefix + "all");
	files.CreateMainFile(subdir, prefix + "other");
	files.CreateMainFile(subdir, prefix + "difference");

	if (isOutputGroups)
	{
		CreateGroupResultFiles(tracks, files, subdir, prefix);
	}
}

void TracerBackScatterPoint::CreateGroupResultFiles(const Tracks &tracks,
													ScatteringFiles &files,
													const string &subdir,
													const string &prefix)
{
	for (size_t i = 0; i < tracks.size(); ++i)
	{
		string groupName = tracks[i].CreateGroupName();
		string filename = prefix + groupName;
		files.CreateGroupFile(subdir, filename);
	}
}

void TracerBackScatterPoint::AllocGroupMatrices(vector<Arr2D> &mtrcs, size_t maxGroupID)
{
	for (size_t i = 0; i < maxGroupID; ++i)
	{
		Arr2D m(1, 1, 4, 4);
		mtrcs.push_back(m);
	}
}
