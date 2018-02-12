#include "TracerBackScatterPoint.h"
#include "global.h"
#include "ScatteringFiles.h"
#include <iostream>

#define BEAM_DIR_LIM		0.9396

using namespace std;

class PointContribution
{
public:
	PointContribution(size_t groupNumber, double normFactor)
		: m_groupNumber(groupNumber),
		  m_normFactor(normFactor)
	{
		groupJones.resize(groupNumber);
		groupMuellers.resize(groupNumber);
		ResetJones();
	}

	void AddToMueller(const Matrix2x2c &jones)
	{
		MuellerMatrix m(jones);
		m *= m_normFactor;
		rest += m;
	}

	void AddToGroup(const Matrix2x2c &jones, size_t groupID)
	{
		groupJones[groupID] += jones;
	}

	void SumGroupTotal()
	{
		for (size_t gr = 0; gr < m_groupNumber; ++gr)
		{
			MuellerMatrix m(groupJones[gr]);
			m *= m_normFactor;
			groupMuellers[gr] += m;
			groupTotal += m;
		}

		ResetJones();
	}

	void SumTotal()
	{
		total += groupTotal;
		total += rest;
	}

	const MuellerMatrix &GetGroupTotal() const
	{
		return groupTotal;
	}

	const MuellerMatrix &GetTotal() const
	{
		return total;
	}

	const MuellerMatrix &GetRest() const
	{
		return rest;
	}

	const MuellerMatrix &GetGroupMueller(size_t groupID)
	{
		return groupMuellers.at(groupID);
	}

	void Reset()
	{
		groupTotal.Reset();
		rest.Reset();
		total.Reset();

		for (MuellerMatrix &m : groupMuellers)
		{
			m.Reset();
		}
	}

private:
	std::vector<Matrix2x2c> groupJones;
	std::vector<MuellerMatrix> groupMuellers;

	MuellerMatrix groupTotal;
	MuellerMatrix rest;
	MuellerMatrix total;

	double m_groupNumber;
	double m_normFactor;

	void ResetJones()
	{
		for (Matrix2x2c &j : groupJones)
		{
			j.Fill(0.f);
		}
	}
};

TracerBackScatterPoint::TracerBackScatterPoint(Particle *particle, int reflNum,
											   const std::string &resultFileName)
	: Tracer(particle, reflNum, resultFileName)
{

}

void TracerBackScatterPoint::Trace(const AngleRange &betaRange, const AngleRange &gammaRange,
								   const Tracks &tracks, double wave)
{
	size_t groupNum = tracks.size();
	normIndex = gammaRange.step/gammaRange.norm;

	m_wavelength = wave;
	CalcTimer timer;
	long long count = 0;

	string dir = CreateFolder(m_resultDirName);
	string tableHead = GetTableHead(betaRange);

	string fulldir = dir + m_resultDirName + '\\';
	ofstream logfile(fulldir + "log.txt", ios::out);

	ScatteringFiles resFiles(dir, m_resultDirName, tableHead);
	PointContribution originContrib(groupNum, normIndex);
	CreateDir(fulldir + "res");
	CreateResultFiles(resFiles, tracks, "res");

	vector<Arr2D> groupResultM;
	AllocGroupMatrices(groupResultM, groupNum);

	ScatteringFiles corFiles(dir, m_resultDirName, tableHead);
	PointContribution correctedContrib(groupNum, normIndex);
	CreateDir(fulldir + "cor");
	CreateResultFiles(corFiles, tracks, "cor", "cor_");

	vector<Arr2D> groupResultM_cor;
	AllocGroupMatrices(groupResultM_cor, groupNum);

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
			m_tracing->SplitBeamByParticle(beta, gamma, outBeams);

			m_incomingEnergy += m_tracing->GetIncomingEnergy();

			HandleBeams(outBeams, tracks, originContrib, correctedContrib);
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

		OutputContribution(groupNum, originContrib, resFiles, degBeta);
		OutputContribution(groupNum, correctedContrib, corFiles, degBeta, "cor_");
	}

	EraseConsoleLine(50);
	cout << "100%";

	if (isOutputGroups)
	{
		for (size_t group = 0; group < groupNum; ++group)
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

void TracerBackScatterPoint::OutputContribution(size_t groupNumber, PointContribution &contrib,
								ScatteringFiles &files, double degree, string prefix)
{
	contrib.SumTotal();

	ofstream *all = files.GetMainFile(prefix + "all");
	*(all) << degree << ' ' << m_incomingEnergy << ' ';
	*(all) << contrib.GetTotal() << endl;

	if (isOutputGroups)
	{
		for (size_t gr = 0; gr < groupNumber; ++gr)
		{
			ofstream &file = *(files.GetGroupFile(gr));
			file << degree << ' ' << m_incomingEnergy << ' ';
			file << contrib.GetGroupMueller(gr) << endl;
		}
	}

	if (isCalcOther)
	{
		ofstream &other = *(files.GetMainFile(prefix + "other"));
		other << degree << ' ' << m_incomingEnergy << ' ';
		other << contrib.GetRest() << endl;

		ofstream &diff = *(files.GetMainFile(prefix + "difference"));
		diff << degree << ' ' << m_incomingEnergy << ' ';
		diff << contrib.GetGroupTotal() << endl;
	}

	contrib.Reset();
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
							const string &subdir, const string &prefix)
{
	files.CreateMainFile(subdir, prefix + "all");
	files.CreateMainFile(subdir, prefix + "other");
	files.CreateMainFile(subdir, prefix + "difference");

	if (isOutputGroups)
	{
		CreateGroupResultFiles(tracks, files, prefix);
	}
}

void TracerBackScatterPoint::CreateGroupResultFiles(const Tracks &tracks, ScatteringFiles &files,
									const string &prefix)
{
	for (size_t i = 0; i < tracks.size(); ++i)
	{
		string groupName = tracks[i].CreateGroupName();
		string filename = prefix + groupName;
		files.CreateGroupFile("", filename);
	}
}

void TracerBackScatterPoint::HandleBeams(std::vector<Beam> &outBeams, const Tracks &tracks,
										 PointContribution &originContrib,
										 PointContribution &correctedContrib)
{
	Point3d vr(0, 0, 1);
	Point3d vf = -m_polarizationBasis;

	for (Beam &beam : outBeams)
	{
		if (beam.direction.cz < BEAM_DIR_LIM)
		{
			continue;
		}

		int groupID = tracks.FindGroup(beam.id);

		if (groupID < 0 && !isCalcOther)
		{
			continue;
		}

		beam.RotateSpherical(-m_incidentDir, m_polarizationBasis);

		Point3f T = CrossProduct(beam.e, beam.direction);
		T = T/Length(T); // basis of beam

		Point3f center = beam.Center();
		double lng_proj0 = beam.opticalPath + DotProduct(center, beam.direction);

		matrixC Jx(2, 2);
		CalcMultiplyOfJmatrix(beam, T, vf, vr, lng_proj0, Jx);

		// correction
		matrixC Jx_cor = Jx;
		Jx_cor[0][1] -= Jx_cor[1][0];
		Jx_cor[0][1] /= 2;
		Jx_cor[1][0] = -Jx_cor[0][1];

		Matrix2x2c _Jx(Jx), _Jx_cor(Jx_cor);

		if (groupID < 0 && isCalcOther)
		{
			originContrib.AddToMueller(_Jx);
			correctedContrib.AddToMueller(_Jx_cor);
		}
		else
		{
			originContrib.AddToGroup(_Jx, groupID);
			correctedContrib.AddToGroup(_Jx_cor, groupID);
		}
	}

	originContrib.SumGroupTotal();
	correctedContrib.SumGroupTotal();
}

void TracerBackScatterPoint::AllocGroupMatrices(vector<Arr2D> &mtrcs, size_t maxGroupID)
{
	for (size_t i = 0; i < maxGroupID; ++i)
	{
		Arr2D m(1, 1, 4, 4);
		mtrcs.push_back(m);
	}
}
