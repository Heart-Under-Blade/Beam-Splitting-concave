#include "Tracer.h"

#include <ostream>
#include <iostream>
#include <iomanip>
#include <map>
#include <assert.h>
#include "global.h"
#include "macro.h"
//#ifdef _TRACK_ALLOW
//std::ofstream trackMapFile("tracks_deb.dat", std::ios::out);
//#endif

#include "TracingConvex.h"
#include "TracingConcave.h"

#define BEAM_DIR_LIM		0.9396
#define SPHERE_RING_NUM		180		// number of rings no the scattering sphere
#define BIN_SIZE			M_PI/SPHERE_RING_NUM
#define RES_EXT ".dat"

using namespace std;

class ContributionFiles : public map<string, ofstream*>
{
public:
	ContributionFiles(const string &dir)
		: map<string, ofstream*>()
	{
		m_dir = dir;
		int i = m_dir.size()-2;

		while (i >= 0 && dir[i] != '\\')
		{
			--i;
		}

		int from = i+1;
		int num = (m_dir.size()-from)-1;
		m_folder = m_dir.substr(from, num);
	}

	void CreateFile(const string &subdir, const string &name)
	{
		string filename = m_dir + subdir + '\\' + name + "__" + m_folder + RES_EXT;
		ofstream *file = new ofstream(filename, ios::out);
		assert(file->is_open());

		(*file) << setprecision(numeric_limits<long double>::digits10 + 1);
		this->insert(pair<string, ofstream*>(name, file));
	}

	~ContributionFiles()
	{
		for (const auto &p : (*this))
		{
			p.second->close();
			delete p.second;
		}
	}

private:
	string m_dir;
	string m_folder;
};

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

Tracer::Tracer(Particle *particle, int reflNum, const string &resultFileName)
	: m_incidentDir(0, 0, -1), // down direction
	  m_polarizationBasis(0, 1, 0),
	  m_resultDirName(resultFileName)
{
	if (particle->IsComplicated())
	{
		m_tracing = new TracingConcave(particle, m_incidentDir, true,
									   m_polarizationBasis, reflNum);
	}
	else
	{
		m_tracing = new TracingConvex(particle, m_incidentDir, true,
									  m_polarizationBasis, reflNum);
	}

	m_symmetry = m_tracing->m_particle->GetSymmetry();
}

Tracer::~Tracer()
{
	delete m_tracing;
}

void OutputOrientationToLog(int i, int j, ostream &logfile)
{
	logfile << "i: " << i << ", j: " << j << endl;
	logfile.flush();
}

void Tracer::WriteResultToSeparateFilesGO(double NRM, int EDF, const string &dir,
										  const Tracks &tracks)
{
	for (size_t i = 0; i < m_sepatateMatrices.size(); ++i)
	{
		if (tracks[i].size != 0)
		{
			string subname = tracks[i].CreateGroupName();
			ExtractPeaksGO(EDF, NRM, m_sepatateMatrices[i]);
			WriteResultsToFileGO(NRM, dir + subname, m_sepatateMatrices[i]);
		}
	}
}

void Tracer::TraceIntervalGO(int betaNumber, int gammaNumber, const Tracks &tracks)
{
	CalcTimer timer;
	string dirName = CreateDir(m_resultDirName);

#ifdef _CHECK_ENERGY_BALANCE
	m_incomingEnergy = 0;
	m_outcomingEnergy = 0;
#endif

	double betaNorm = m_symmetry.beta/betaNumber;
	double gammaNorm = m_symmetry.gamma/gammaNumber;

	m_sepatateMatrices.resize(tracks.size());

//	m_totalMtrx.scatMatrix = Arr2D(1, thetaNum+1, 4, 4);
//	m_totalMtrx.scatMatrix.ClearArr();

	vector<Beam> outBeams;
	double beta, gamma;

	m_startTime = timer.Start();
	cout << "Started at " << ctime(&m_startTime) << endl;

	for (int i = 0; i < betaNumber; ++i)
	{
		beta = (i + 0.5)*betaNorm;

		for (int j = 0; j < gammaNumber; ++j)
		{
			gamma = (j + 0.5)*gammaNorm;
			m_tracing->SplitBeamByParticle(beta, gamma, outBeams);

#ifdef _CHECK_ENERGY_BALANCE
//			m_incomingEnergy += m_tracing->GetIncomingEnergy()*sin(beta);
#endif
			HandleBeamsGO(outBeams, beta, tracks);
			outBeams.clear();
//			OutputOrientationToLog(i, j, logfile);
		}

		OutputProgress(betaNumber, i, timer);
	}

	double D_tot = CalcTotalScatteringEnergy();
	long long orNum = gammaNumber * betaNumber;
	double NRM = CalcNorm(orNum);

#ifdef _CHECK_ENERGY_BALANCE
//	for (int deg = 0; deg <= thetaNum; ++deg)
//	{
//		m_outcomingEnergy += m_muller(0, deg, 0, 0)*NRM;
//	}
#endif

	int EDF = 0;
	WriteResultToSeparateFilesGO(NRM, EDF, dirName, tracks);

	ExtractPeaksGO(EDF, NRM, m_totalMtrx);
	WriteResultsToFileGO(NRM, dirName + "all", m_totalMtrx);
	OutputStatisticsGO(orNum, D_tot, NRM, timer);
}

void Tracer::OutputProgress(int betaNumber, long long count, CalcTimer &timer)
{
	EraseConsoleLine(50);
	cout << (count*100)/(betaNumber+1) << '%'
		 << '\t' << timer.Elapsed();
}

void Tracer::ExtractPeaksGO(int EDF, double NRM, ContributionGO &contr)
{
	//Analytical averaging over alpha angle
	double b[3], f[3];
	b[0] =  contr.back[0][0];
	b[1] = (contr.back[1][1] - contr.back[2][2])/2.0;
	b[2] =  contr.back[3][3];

	f[0] =  contr.forw[0][0];
	f[1] = (contr.forw[1][1] + contr.forw[2][2])/2.0;
	f[2] =  contr.forw[3][3];

	// Extracting the forward and backward peak in a separate file if needed
	if (EDF)
	{
		std::ofstream bck("back.dat", std::ios::out);
		std::ofstream frw("forward.dat", std::ios::out);
		frw << "M11 M22/M11 M33/M11 M44/M11";
		bck << "M11 M22/M11 M33/M11 M44/M11";

		if (f[0] <= DBL_EPSILON)
		{
			frw << "\n0 0 0 0";
		}
		else
		{
			frw << "\n" << f[0]*NRM
				<< " " << f[1]/f[0]
				<< " " << f[1]/f[0]
				<< " " << f[2]/f[0];
		}

		if (b[0] <= DBL_EPSILON)
		{
			bck << "\n0 0 0 0";
		}
		else
		{
			bck << "\n" << b[0]*NRM
				<< " " << b[1]/b[0]
				<< " " << -b[1]/b[0]
				<< " " << b[2]/b[0];
		}

		bck.close();
		frw.close();
	}
	else
	{
		contr.scatMatrix(0,SPHERE_RING_NUM,0,0) += f[0];
		contr.scatMatrix(0,0,0,0) += b[0];
		contr.scatMatrix(0,SPHERE_RING_NUM,1,1) += f[1];
		contr.scatMatrix(0,0,1,1) += b[1];
		contr.scatMatrix(0,SPHERE_RING_NUM,2,2) += f[1];
		contr.scatMatrix(0,0,2,2) -= b[1];
		contr.scatMatrix(0,SPHERE_RING_NUM,3,3) += f[2];
		contr.scatMatrix(0,0,3,3) += b[2];
	}
}

void Tracer::WriteResultsToFileGO(double NRM, const string &filename,
								  ContributionGO &contr)
{
	string name = GetUniqueFileName(filename);
	ofstream allFile(name, std::ios::out);

	allFile << "tetta M11 M12/M11 M21/M11 M22/M11 M33/M11 M34/M11 M43/M11 M44/M11";

	for (int j = SPHERE_RING_NUM; j >= 0; j--)
	{
		double tmp0 = 180.0/SPHERE_RING_NUM*(SPHERE_RING_NUM-j);
		double tmp1 = (j == 0) ? -(0.25*180.0)/SPHERE_RING_NUM : 0;
		double tmp2 = (j == (int)SPHERE_RING_NUM) ? (0.25*180.0)/SPHERE_RING_NUM : 0;

		// Special case in first and last step
		allFile << '\n' << tmp0 + tmp1 + tmp2;

		double sn = (j == 0 || j == (int)SPHERE_RING_NUM)
				? 1-cos(BIN_SIZE/2.0)
				: (cos((j-0.5)*BIN_SIZE)-cos((j+0.5)*BIN_SIZE));

		matrix bf = contr.scatMatrix(0, j);

		if (bf[0][0] <= DBL_EPSILON)
		{
			allFile << " 0 0 0 0 0 0 0 0";
		}
		else
		{
			allFile << ' ' << bf[0][0]*NRM/(2.0*M_PI*sn)
					<< ' ' << bf[0][1]/bf[0][0]
					<< ' ' << bf[1][0]/bf[0][0]
					<< ' ' << bf[1][1]/bf[0][0]
					<< ' ' << bf[2][2]/bf[0][0]
					<< ' ' << bf[2][3]/bf[0][0]
					<< ' ' << bf[3][2]/bf[0][0]
					<< ' ' << bf[3][3]/bf[0][0];
		}
	}

	allFile.close();
}

void Tracer::TraceRandomPO(int betaNumber, int gammaNumber, const Cone &bsCone,
						   const Tracks &tracks, double wave)
{
	m_wavelength = wave;
	CalcTimer timer;
	long long count = 0;

	Arr2D M(bsCone.phiCount+1, bsCone.thetaCount+1, 4, 4);
	ofstream outFile(m_resultDirName, ios::out);

	vector<Beam> outBeams;
	double beta, gamma;

	double betaStep = m_symmetry.beta/betaNumber;
	double gammaStep = m_symmetry.gamma/gammaNumber;

	int halfGammaCount = gammaNumber/2;
	gNorm = gammaStep/m_symmetry.gamma;

	timer.Start();

	for (int i = 0; i <= betaNumber; ++i)
	{
		beta = i*betaStep;

		for (int j = -halfGammaCount; j <= halfGammaCount; ++j)
		{
			gamma = j*gammaStep;

			m_tracing->SplitBeamByParticle(beta, gamma, outBeams);

			CleanJ(tracks.size(), bsCone);
			HandleBeamsPO(outBeams, bsCone, tracks);

			outBeams.clear();
			AddResultToMatrix(M, bsCone, gNorm);
		}

		WriteConusMatrices(outFile, M, bsCone);

		OutputProgress(betaNumber, count, timer);
		++count;
	}

	outFile.close();
}

void Tracer::AllocGroupMatrices(vector<Arr2D> &mtrcs, size_t maxGroupID)
{
	for (size_t i = 0; i < maxGroupID; ++i)
	{
		Arr2D m(1, 1, 4, 4);
		mtrcs.push_back(m);
	}
}

void Tracer::CreateGroupResultFiles(const string &tableHead, const Tracks &tracks,
									const string &dirName,
									vector<ofstream*> &groupFiles,
									const string &prefix)
{
	for (size_t i = 0; i < tracks.size(); ++i)
	{
		string groupName = tracks[i].CreateGroupName();
		string filename = dirName + prefix + groupName + "__" + m_resultDirName + ".dat";
		ofstream *file = new ofstream(filename, ios::out);
		(*file) << tableHead;
		(*file) << std::setprecision(std::numeric_limits<long double>::digits10 + 1);
		groupFiles.push_back(file);
	}
}

void Tracer::AllocJ(std::vector<Arr2DC> &j, int m, int n, int size)
{
	j.clear();
	Arr2DC tmp(m, n, 2, 2);
	tmp.ClearArr();

	for(int i = 0; i < size; i++)
	{
		j.push_back(tmp);
	}
}

void Tracer::CleanJ(std::vector<Arr2DC> &j)
{
	for (Arr2DC &m : j)
	{
		m.ClearArr();
	}
}

void Tracer::OutputToGroupFiles(double degBeta, std::vector<ofstream *> &groupFiles,
								std::vector<Arr2D> &groupResultM, size_t size)
{
	for (size_t group = 0; group < size; ++group)
	{
		Arr2D &mtrx = groupResultM[group];
		matrix m1 = mtrx(0, 0);
		ofstream &file = *(groupFiles[group]);
		file << degBeta << ' ' << m_incomingEnergy << ' ';
		file << m1 << endl;
	}

	groupResultM.clear();
	AllocGroupMatrices(groupResultM, size);
}

void Tracer::OutputToAllFile(ofstream &diffFile, ofstream &otherFile,
							 double degBeta, ofstream &allFile, Arr2D &all,
							 Arr2D &other)
{
	allFile << degBeta << ' ' << m_incomingEnergy << ' ';
	matrix m = all(0, 0);
	allFile << m << endl;
	all.ClearArr();

	if (isCalcOther)
	{
		otherFile << degBeta << ' ' << m_incomingEnergy << ' ';
		matrix m0 = other(0, 0);
		otherFile << m0 << endl;
		other.ClearArr();

		diffFile << degBeta << ' ' << m_incomingEnergy << ' ';
		matrix m00 = m - m0;
		diffFile << m00 << endl;
	}
}

void Tracer::OutputStatisticsPO(CalcTimer &timer, long long orNumber, const string &path)
{
	string startTime = ctime(&m_startTime);
	string totalTime = timer.Elapsed();
	time_t end = timer.Stop();
	string endTime = ctime(&end);

	m_statistics += "\nStart of calculation = " + startTime
			+ "End of calculation   = " + endTime
			+ "\nTotal time of calculation = " + totalTime
			+ "\nTotal number of body orientation = " + to_string(orNumber);

	if (isNanOccured)
	{
		m_statistics += "\n\nWARNING! NAN values occured. See 'log.txt'";
	}

	ofstream out(path + "\\out.dat", ios::out);

	out << m_statistics;
	out.close();

	cout << m_statistics;
}

void Tracer::setIsOutputGroups(bool value)
{
	isOutputGroups = value;
}

void Tracer::OutputContribution(size_t groupNumber, PointContribution &contrib,
								vector<ofstream*> groupFiles, ContributionFiles &files,
								double degree, string prefix)
{
	contrib.SumTotal();

	ofstream &all = *(files[prefix + "all"]);
	all << degree << ' ' << m_incomingEnergy << ' ';
	all << contrib.GetTotal() << endl;

	if (isOutputGroups)
	{
		for (size_t gr = 0; gr < groupNumber; ++gr)
		{
			ofstream &file = *(groupFiles[gr]);
			file << degree << ' ' << m_incomingEnergy << ' ';
			file << contrib.GetGroupMueller(gr) << endl;
		}
	}

	if (isCalcOther)
	{
		ofstream &other = *(files[prefix + "other"]);
		other << degree << ' ' << m_incomingEnergy << ' ';
		other << contrib.GetRest() << endl;

		ofstream &diff = *(files[prefix + "difference"]);
		diff << degree << ' ' << m_incomingEnergy << ' ';
		diff << contrib.GetGroupTotal() << endl;
	}

	contrib.Reset();
}

string Tracer::GetTableHead(const AngleRange &range)
{
	return to_string(range.number) + ' '
			+ to_string(RadToDeg(range.max)) + ' '
			+ to_string(RadToDeg(range.step)) + '\n'
			+ "beta cr_sec M11 M12 M13 M14 M21 M22 M23 M24 M31 M32 M33 M34 M41 M42 M43 M44"
			+ '\n';
}

void Tracer::TraceBackScatterPointPO(const AngleRange &betaRange, const AngleRange &gammaRange,
									 const Tracks &tracks, double wave)
{
	int groupNumber = tracks.size();
	gNorm = gammaRange.step/gammaRange.norm;

	PointContribution originContrib(groupNumber, gNorm);
	PointContribution correctedContrib(groupNumber, gNorm);

	m_wavelength = wave;
	CalcTimer timer;
	long long count = 0;

	string dir = CreateDir(m_resultDirName);

	ofstream logfile(dir + "\\log.txt", ios::out);
	vector<ofstream*> groupFiles;

	string resDir = CreateDir2(dir + "res");
	string corDir = CreateDir2(dir + "cor");

	string tableHead = GetTableHead(betaRange);

	if (isOutputGroups)
	{
		CreateGroupResultFiles(tableHead, tracks, resDir, groupFiles);
	}

	ContributionFiles files(dir);
	files.CreateFile("res", "all");
	files.CreateFile("res", "other");
	files.CreateFile("res", "difference");
	files.CreateFile("cor", "cor_all");
	files.CreateFile("cor", "cor_other");
	files.CreateFile("cor", "cor_difference");

	for (const auto &p : files)
	{
		*(p.second) << tableHead;
	}

	vector<ofstream*> groupFiles_cor;

	if (isOutputGroups)
	{
		CreateGroupResultFiles(tableHead, tracks, corDir, groupFiles_cor, "cor_");
	}

	vector<Arr2D> groupResultM, groupResultM_cor;
	AllocGroupMatrices(groupResultM, groupNumber);
	AllocGroupMatrices(groupResultM_cor, groupNumber);

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

			HandleBeamsBackScatterPO(outBeams, tracks, originContrib, correctedContrib);
			outBeams.clear();
			OutputOrientationToLog(i, j, logfile);

			if (isNan)
			{
				logfile << "------------------";
				isNan = false;
				continue;
			}
		}

		m_incomingEnergy *= gNorm;

		double degBeta = RadToDeg(beta);

		OutputContribution(groupNumber, originContrib, groupFiles, files,
						   degBeta);

		OutputContribution(groupNumber, correctedContrib, groupFiles_cor, files,
						   degBeta, "cor_");
	}

	EraseConsoleLine(50);
	cout << "100%";

	if (isOutputGroups)
	{
		for (size_t group = 0; group < groupResultM.size(); ++group)
		{
			ofstream &file = *(groupFiles[group]);
			file.close();

			ofstream &file_cor = *(groupFiles_cor[group]);
			file_cor.close();
		}
	}

	long long orNumber = betaRange.number * gammaRange.number;
	OutputStatisticsPO(timer, orNumber, m_resultDirName);
}

//REF: объединить с предыдущим
void Tracer::TraceIntervalPO2(int betaNumber, int gammaNumber, const Cone &bsCone,
							  const Tracks &tracks, double wave)
{
	m_wavelength = wave;
	CalcTimer timer;
	long long count = 0;

	Arr2D M(bsCone.phiCount+1, bsCone.thetaCount+1, 4, 4);
	ofstream outFile(m_resultDirName, ios::out);

	vector<Beam> outBeams;
	double beta, gamma;

	double betaNorm = m_symmetry.beta/betaNumber;
	double gammaNorm = m_symmetry.gamma/gammaNumber;

	int halfGammaCount = gammaNumber/2;
	gNorm = gammaNorm/m_symmetry.gamma;

	timer.Start();

	for (int i = 0; i <= betaNumber; ++i)
	{
		beta = i*betaNorm;

		for (int j = -halfGammaCount; j <= halfGammaCount; ++j)
		{
			gamma = j*gammaNorm;

			for (size_t groupID = 0; groupID < tracks.size(); ++groupID)
			{
				m_tracing->SplitBeamByParticle(beta, gamma, tracks[groupID].tracks, outBeams);

				CleanJ(tracks.size(), bsCone);
				HandleBeamsPO2(outBeams, bsCone, groupID);

				outBeams.clear();
				AddResultToMatrix(M, bsCone, gNorm);
			}
		}

		WriteConusMatrices(outFile, M, bsCone);

		OutputProgress(betaNumber, count, timer);
		++count;
	}

	outFile.close();
}

void Tracer::WriteConusMatrices(ofstream &outFile, const Arr2D &sum,
								const Cone &bsCone)
{
	outFile << to_string(bsCone.radius)	<< ' '
			<< to_string(bsCone.thetaCount) << ' '
			<< to_string(bsCone.phiCount+1);

	for (int t = 0; t <= bsCone.thetaCount; ++t)
	{
		double tt = RadToDeg(t*bsCone.dTheta);

		for (int p = 0; p <= bsCone.phiCount; ++p)
		{
			double fi = -((double)p)*bsCone.dPhi;
			double degPhi = RadToDeg(-fi);
			outFile << endl << tt << " " << degPhi << " ";

			matrix m = sum(p, t);
			outFile << m;
		}
	}
}

void Tracer::AddResultToMatrix(Arr2D &M, const Cone &bsCone, double norm)
{
	for (size_t q = 0; q < J.size(); ++q)
	{
		for (int t = 0; t <= bsCone.thetaCount; ++t)
		{
			for (int p = 0; p <= bsCone.phiCount; ++p)
			{
				matrix m = Mueller(J[q](p, t));
				m *= norm;
				M.insert(p, t, m);
			}
		}
	}
}

void Tracer::AddResultToMatrix(Arr2D &M, std::vector<Arr2DC> &j, double norm)
{
	for (size_t q = 0; q < j.size(); ++q)
	{
		matrix m = Mueller(j[q](0, 0));
		m *= norm;
		M.insert(0, 0, m);
	}
}

void Tracer::AddResultToMatrices(vector<Arr2D> &M, const Cone &bsCone,
								 double norm)
{
	for (size_t q = 0; q < J.size(); ++q)
	{
		for (int t = 0; t <= bsCone.thetaCount; ++t)
		{
			for (int p = 0; p <= bsCone.phiCount; ++p)
			{
				matrix m = Mueller(J[q](0, 0));
				m *= norm;
				M[q].insert(p, t, m);
			}
		}
	}
}

void Tracer::AddResultToMatrices(std::vector<Arr2D> &M)
{
	for (size_t q = 0; q < J.size(); ++q)
	{
		matrix m = Mueller(J[q](0, 0));
		m *= gNorm;
		M[q].insert(0, 0, m);
	}
}
// REF: merge with previous
void Tracer::AddResultToMatrices_cor(std::vector<Arr2D> &M, double norm)
{
	for (size_t q = 0; q < J_cor.size(); ++q)
	{
		matrix m = Mueller(J_cor[q](0, 0));
		m *= norm;
		M[q].insert(0, 0, m);
	}
}

void Tracer::CleanJ(int size, const Cone &bsCone)
{
	J.clear();
	Arr2DC tmp(bsCone.phiCount+1, bsCone.thetaCount+1, 2, 2);
	tmp.ClearArr();

	for(int q = 0; q < size; q++)
	{
		J.push_back(tmp);
	}
}

void Tracer::HandleBeamsGO(std::vector<Beam> &outBeams, double beta)
{
	double sinBeta = sin(beta);

	for (Beam &beam : outBeams)
	{
		beam.RotateSpherical(-m_incidentDir, m_polarizationBasis);

		double cross = m_tracing->BeamCrossSection(beam);
		double area = cross*sinBeta;
		matrix bf = Mueller(beam.J);

		AddToResultMullerGO(beam.direction, bf, area, m_totalMtrx);
	}
}

void Tracer::RotateMuller(const Point3f &dir, matrix &bf)
{
	const float &x = dir.cx;
	const float &y = dir.cy;

	double tmp = y*y;

	tmp = acos(x/sqrt(x*x+tmp));

	if (y < 0)
	{
		tmp = M_2PI-tmp;
	}

	tmp *= -2.0;
	RightRotateMueller(bf, cos(tmp), sin(tmp));
}

void Tracer::AddToResultMullerGO(const Point3f &dir, matrix &bf, double area,
								 ContributionGO &contr)
{
	const float &z = dir.cz;

	// Collect the beam in array
	if (z >= 1-DBL_EPSILON)
	{
		contr.back += area*bf;
	}
	else if (z <= DBL_EPSILON-1)
	{
		contr.forw += area*bf;
	}
	else
	{
		const float &y = dir.cy;

		if (y*y > DBL_EPSILON)
		{	// rotate the Mueller matrix of the beam to appropriate coordinate system
			RotateMuller(dir, bf);
		}

#ifdef _CALC_AREA_CONTIBUTION_ONLY
		bf = matrix(4,4);
		bf.Identity();
#endif
		const unsigned int zenAng = round(acos(z)/(M_PI/SPHERE_RING_NUM));
		contr.scatMatrix.insert(0, zenAng, area*bf);
	}
}

void Tracer::HandleBeamsGO(std::vector<Beam> &outBeams, double beta,
						   const Tracks &tracks)
{
	double sinBeta = sin(beta);

	for (Beam &beam : outBeams)
	{
		int grID = tracks.FindGroup(beam.id);

		if (grID < 0)
		{
			continue;
		}

		beam.RotateSpherical(-m_incidentDir, m_polarizationBasis);

		double cross = m_tracing->BeamCrossSection(beam);
		double area = cross*sinBeta;
		matrix bf = Mueller(beam.J);

		// individual contribution
		AddToResultMullerGO(beam.direction, bf, area, m_sepatateMatrices[grID]);

		// total contribution
		AddToResultMullerGO(beam.direction, bf, area, m_totalMtrx);
	}

	outBeams.clear();
}

double Tracer::CalcNorm(long long orNum)
{
	double &symBeta = m_symmetry.beta;
	double tmp = (/*isRandom*/true) ? symBeta : 1.0;
	double dBeta = -(cos(symBeta) - cos(0));
	return tmp/(orNum*dBeta);
}

double Tracer::CalcTotalScatteringEnergy()
{
	double D_tot = m_totalMtrx.back[0][0] + m_totalMtrx.forw[0][0];

	for (int i = 0; i <= SPHERE_RING_NUM; ++i)
	{
		D_tot += m_totalMtrx.scatMatrix(0, i, 0, 0);
	}

	return D_tot;
}

void Tracer::OutputStartTime(CalcTimer &timer)
{
	m_startTime = timer.Start();
	cout << "Started at " << ctime(&m_startTime) << endl;
}

void Tracer::OutputStatisticsGO(int orNumber, double D_tot, double NRM,
								CalcTimer &timer)
{
	string startTime = ctime(&m_startTime);
	string totalTime = timer.Elapsed();
	time_t end = timer.Stop();
	string endTime = ctime(&end);

	m_statistics += "\nStart of calculation = " + startTime
			+ "End of calculation   = " + endTime
			+ "\nTotal time of calculation = " + totalTime
			+ "\nTotal number of body orientation = " + to_string(orNumber)
			+ "\nTotal scattering energy = " + to_string(D_tot);

#ifdef _CHECK_ENERGY_BALANCE
	double normEnergy = m_incomingEnergy * NRM;
	double passedEnergy = (m_outcomingEnergy/normEnergy)*100;

	m_statistics += "\nTotal incoming energy = " + to_string(normEnergy)
			+ "\nTotal outcoming energy = " + to_string(m_outcomingEnergy)
			+ "\nEnergy passed = " + to_string(passedEnergy) + '%';
#endif

	//	out << "\nAveraged cross section = " << incomingEnergy*NRM;
	ofstream out(m_resultDirName+"_out.dat", ios::out);
	out << m_statistics;
	out.close();

	cout << m_statistics;
}

void Tracer::TraceIntervalGO(int betaNumber, int gammaNumber)
{
	int EDF = 0;
	CalcTimer timer;

#ifdef _CHECK_ENERGY_BALANCE
	m_incomingEnergy = 0;
	m_outcomingEnergy = 0;
#endif

//	double dd = M_PI/SPHERE_RING_NUM;

	double betaNorm = m_symmetry.beta/betaNumber;
	double gammaNorm = m_symmetry.gamma/gammaNumber;

	m_totalMtrx.back.Fill(0);
	m_totalMtrx.forw.Fill(0);

	m_totalMtrx.scatMatrix = Arr2D(1, SPHERE_RING_NUM+1, 4, 4);
	m_totalMtrx.scatMatrix.ClearArr();

	vector<Beam> outBeams;
	double beta, gamma;

	OutputStartTime(timer);

	for (int i = 0; i < betaNumber; ++i)
	{
		beta = (i + 0.5)*betaNorm;

		for (int j = 0; j < gammaNumber; ++j)
		{
			gamma = (j + 0.5)*gammaNorm;
			m_tracing->SplitBeamByParticle(beta, gamma, outBeams);

#ifdef _CHECK_ENERGY_BALANCE
			m_incomingEnergy += m_tracing->GetIncomingEnergy()*sin(beta);
#endif
			HandleBeamsGO(outBeams, beta);
			outBeams.clear();
//			OutputOrientationToLog(i, j, logfile);
		}

		OutputProgress(betaNumber, i, timer);
	}

	double D_tot = CalcTotalScatteringEnergy();
	long long orNum = gammaNumber * betaNumber;
	double NRM = CalcNorm(orNum);

#ifdef _CHECK_ENERGY_BALANCE
	for (int deg = 0; deg <= SPHERE_RING_NUM; ++deg)
	{
		m_outcomingEnergy += m_totalMtrx.scatMatrix(0, deg, 0, 0)*NRM;
	}
#endif

	ExtractPeaksGO(EDF, NRM, m_totalMtrx);
	WriteResultsToFileGO(NRM, m_resultDirName + "_all", m_totalMtrx);
	OutputStatisticsGO(orNum, D_tot, NRM, timer);
}

void Tracer::TraceSingleOrGO(const double &beta, const double &gamma,
							 const Tracks &tracks)
{
	int EDF = 0;
	vector<Beam> outBeams;

	m_totalMtrx.back.Fill(0);
	m_totalMtrx.forw.Fill(0);

	m_totalMtrx.scatMatrix = Arr2D(1, SPHERE_RING_NUM+1, 4, 4);
	m_totalMtrx.scatMatrix.ClearArr();

	double b = DegToRad(beta);
	double g = DegToRad(gamma);
	m_tracing->SplitBeamByParticle(b, g, outBeams);

	HandleBeamsGO(outBeams, beta, tracks);

//	double D_tot = CalcTotalScatteringEnergy();

	ExtractPeaksGO(EDF, 1, m_totalMtrx);
	WriteResultsToFileGO(1, m_resultDirName, m_totalMtrx);
//	WriteStatisticsToFileGO(1, D_tot, 1, timer); // TODO: раскомментить
}

void Tracer::TraceSingleOrPO(const double &beta, const double &gamma,
							 const Cone &bsCone, const Tracks &tracks, double wave)
{
	m_wavelength = wave;
	Arr2D M(bsCone.phiCount+1, bsCone.thetaCount+1, 4, 4);
	ofstream outFile(m_resultDirName, ios::out);
	vector<Beam> outBeams;

	double b = DegToRad(beta);
	double g = DegToRad(gamma);
	m_tracing->SplitBeamByParticle(b, g, outBeams);

	CleanJ(tracks.size(), bsCone);
	HandleBeamsPO(outBeams, bsCone, tracks);

	outBeams.clear();
	AddResultToMatrix(M, bsCone);
	WriteConusMatrices(outFile, M, bsCone);
}

void Tracer::setIsCalcOther(bool value)
{
	isCalcOther = value;
}

void Tracer::CalcJnRot(const Beam &beam, const Point3f &T,
					  const Point3d &vf, const Point3d &vr, matrixC &Jn_rot)
{
	Point3f normal = beam.Normal();

	Point3d vt = CrossProductD(vf, vr);
	vt = vt/LengthD(vt);

	Point3f NT = CrossProduct(normal, T);
	Point3f NE = CrossProduct(normal, beam.e);

	Point3d NTd = Point3d(NT.cx, NT.cy, NT.cz);
	Point3d NEd = Point3d(NE.cx, NE.cy, NE.cz);

	Jn_rot[0][0] = -DotProductD(NTd, vf);
	Jn_rot[0][1] = -DotProductD(NEd, vf);
	Jn_rot[1][0] =  DotProductD(NTd, vt);
	Jn_rot[1][1] =  DotProductD(NEd, vt);
}

void Tracer::HandleBeamsPO(vector<Beam> &outBeams, const Cone &bsCone,
						   const Tracks &tracks)
{
	for (Beam &beam : outBeams)
	{
		int groupID = tracks.FindGroup(beam.id);

		if (groupID < 0)
		{
			continue;
		}

		beam.RotateSpherical(-m_incidentDir, m_polarizationBasis);

		Point3f center = beam.Center();
		double lng_proj0 = beam.opticalPath + DotProduct(center, beam.direction);

		Point3f T = CrossProduct(beam.e, beam.direction);
		T = T/Length(T); // basis of beam

		for (int i = 0; i <= bsCone.phiCount; ++i)
		{
			for (int j = 0; j <= bsCone.thetaCount; ++j)
			{	//
				double f = i*bsCone.dPhi;
				double t = j*bsCone.dTheta;

				double sinT = sin(t), sinF = sin(f), cosF = cos(f);

				Point3d vr(sinT*cosF, sinT*sinF, cos(t));
				Point3d vf = (j == 0) ? -m_polarizationBasis
									  : Point3d(-sinF ,cosF ,0);
				// OPT: вышеописанные параметры можно вычислить один раз и занести в массив

				matrixC Jx(0, 0);
				CalcMultiplyOfJmatrix(beam, T, vf, vr, lng_proj0, Jx);
				J[groupID].insert(i, j, Jx);
			}
		}
	}
}

void Tracer::HandleBeamsPO2(vector<Beam> &outBeams, const Cone &bsCone, int groupID)
{
	for (unsigned int i = 0; i < outBeams.size(); ++i)
	{
		Beam &beam = outBeams.at(i);

		beam.RotateSpherical(-m_incidentDir, m_polarizationBasis);

		Point3f center = beam.Center();
		Point3d center_d = Point3d(center);
		double lng_proj0 = beam.opticalPath + DotProduct(center, beam.direction);

		Point3f T = CrossProduct(beam.e, beam.direction);
		T = T/Length(T); // basis of beam

		for (int i = 0; i <= bsCone.phiCount; ++i)
		{
			double f = i*bsCone.dPhi;
			double sinF = sin(f), cosF = cos(f);

			for (int j = 0; j <= bsCone.thetaCount; ++j)
			{
				double t = j*bsCone.dTheta;
				double sinT = sin(t);

				Point3d vr(sinT*cosF, sinT*sinF, cos(t));
				Point3d vf = (j == 0) ? -m_polarizationBasis
									  : Point3d(-sinF, cosF ,0);
				matrixC Jn_rot(2, 2);
				CalcJnRot(beam, T, vf, vr, Jn_rot);

				complex fn(0, 0);
				fn = beam.DiffractionIncline(vr, m_wavelength);

				double dp = DotProductD(vr, center_d);
				complex tmp = exp_im(M_2PI*(lng_proj0-dp)/m_wavelength);
				matrixC fn_jn = beam.J * tmp;

				matrixC c = fn*Jn_rot*fn_jn;
				J[groupID].insert(i, j, c);
			}
		}
	}
}


void Tracer::HandleBeamsBackScatterPO(std::vector<Beam> &outBeams, const Tracks &tracks,
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

void Tracer::CalcMultiplyOfJmatrix(const Beam &beam, const Point3f &T,
								   const Point3d &vf, const Point3d &vr,
								   double lng_proj0, matrixC &Jx)
{
	matrixC Jn_rot(2, 2);
	CalcJnRot(beam, T, vf, vr, Jn_rot);

	complex fn(0, 0);
	fn = beam.DiffractionIncline(vr, m_wavelength);

	if (isnan(real(fn)))
	{
		isNanOccured = isNan = true;
		return;
	}

	double dp = DotProductD(vr, Point3d(beam.Center()));
	complex tmp = exp_im(M_2PI*(lng_proj0-dp)/m_wavelength);
	matrixC fn_jn = beam.J * tmp;

	Jx = fn*Jn_rot*fn_jn;
}
