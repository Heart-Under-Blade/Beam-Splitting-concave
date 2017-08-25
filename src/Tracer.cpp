#include "Tracer.h"

#include <iostream>
#include <assert.h>
#include "global.h"
#include "macro.h"

#define BEAM_DIR_LIM	0.9396
#define SPHERE_RING_NUM 180
#define BIN_SIZE		M_PI/SPHERE_RING_NUM

using namespace std;

Tracer::Tracer(Tracing *tracing, const string resultFileName)
	: m_tracing(tracing),
	  m_incidentDir(0, 0, -1), // down direction
	  m_polarizationBasis(0, 1, 0),
	  m_resultDirName(resultFileName)
{
	m_symmetry = m_tracing->m_particle->GetSymmetry();
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
	int EDF = 0;
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
			OutputState(i, j);
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

	WriteResultToSeparateFilesGO(NRM, EDF, dirName, tracks);

	ExtractPeaksGO(EDF, NRM, m_totalMtrx);
	WriteResultsToFileGO(NRM, dirName + "all", m_totalMtrx);

	OutputStatisticsGO(orNum, D_tot, NRM, timer);
}

void Tracer::OutputProgress(int betaNumber, long long count, CalcTimer &timer)
{
	EraseConsoleLine(50);
	cout << (count*100)/betaNumber << '%'
		 << '\t' << timer.Elapsed();
}

void Tracer::ExtractPeaksGO(int EDF, double NRM, Contribution &contr)
{
	//Analytical averaging over alpha angle
	double b[3], f[3];
	b[0] = contr.back[0][0];
	b[1] = (contr.back[1][1] - contr.back[2][2])/2.0;
	b[2] = contr.back[3][3];

	f[0] = contr.forw[0][0];
	f[1] = (contr.forw[1][1] + contr.forw[2][2])/2.0;
	f[2] = contr.forw[3][3];

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
								  Contribution &contr)
{
	string name = GetDATFileName(filename);
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

string Tracer::GetDATFileName(const string &filename)
{
	string name = filename + ".dat";

	for (int i = 1; ifstream(name) != NULL; ++i)
	{
		name = filename + '(' + to_string(i) + ')' + ".dat";
	}

	return name;
}

void Tracer::TraceRandomPO(int betaNumber, int gammaNumber, const Cone &bsCone,
						   const Tracks &tracks, double wave)
{
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
			HandleBeamsPO(outBeams, bsCone, wave, tracks);

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

void Tracer::CreateGroupResultFiles(const AngleRange &betaRange, const Tracks &tracks,
									const string &dirName,
									vector<ofstream*> &groupFiles)
{
	for (size_t i = 0; i < tracks.size(); ++i)
	{
		string groupName = tracks[i].CreateGroupName();
		string filename = dirName + groupName + "__" + m_resultDirName + ".dat";
		ofstream *file = new ofstream(filename, ios::out);
		OutputTableHead(betaRange, *file);
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

void Tracer::OutputTableHead(const AngleRange &betaRange, ofstream &allFile)
{
	allFile << betaRange.number << ' '
			<< RadToDeg(betaRange.max) << ' '
			<< RadToDeg(betaRange.step) << endl;

	allFile << "beta cr_sec M11 M12 M13 M14 M21 M22 M23 M24 M31 M32 M33 M34 M41 M42 M43 M44"
			<< endl;
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

void Tracer::CreateResultFile(ofstream &file, const string &dirName, const string &fileName,
							  const AngleRange &betaRange)
{
	file.open(dirName + fileName + "__" + m_resultDirName + ".dat", ios::out);
	OutputTableHead(betaRange, file);
}

void Tracer::CreateResultFiles(ofstream &all, ofstream &diff, ofstream &other,
							   const AngleRange &betaRange, string dirName, Arr2D &otherArr)
{
	CreateResultFile(all, dirName, "all", betaRange);

	if (isCalcOther)
	{
		CreateResultFile(other, dirName, "other", betaRange);
		CreateResultFile(diff, dirName, "difference", betaRange);
		otherArr = Arr2D(1, 1, 4, 4);
	}
}

void Tracer::OutputStatisticsPO(CalcTimer &timer, long long orNumber)
{
	string startTime = ctime(&m_startTime);
	string totalTime = timer.Elapsed();
	time_t end = timer.Stop();
	string endTime = ctime(&end);

	m_statistics += "\nStart of calculation = " + startTime
			+ "End of calculation   = " + endTime
			+ "\nTotal time of calculation = " + totalTime
			+ "\nTotal number of body orientation = " + to_string(orNumber);

	ofstream out("out.dat", ios::out);
	out << m_statistics;
	out.close();

	cout << m_statistics;
}

void Tracer::TraceBackScatterPointPO(const AngleRange &betaRange, const AngleRange &gammaRange,
									 const Tracks &tracks, double wave)
{
	CalcTimer timer;
	long long count = 0;

	/*string dirName = */CreateDir(m_resultDirName);

	ofstream allFile, otherFile, diffFile;
	vector<ofstream*> groupFiles;

	string resDirName = CreateDir(m_resultDirName + "\\res");
	CreateGroupResultFiles(betaRange, tracks, resDirName, groupFiles);
	CreateResultFiles(allFile, diffFile, otherFile, betaRange, resDirName, Other);

	ofstream allFile_cor, otherFile_cor, diffFile_cor;
	vector<ofstream*> groupFiles_cor;

	string corDirName = CreateDir(m_resultDirName + "\\cor");
	CreateGroupResultFiles(betaRange, tracks, corDirName, groupFiles_cor);
	CreateResultFiles(allFile_cor, diffFile_cor, otherFile_cor, betaRange, corDirName, Other_cor);

	AllocJ(J, 1, 1, tracks.size());
	AllocJ(J_cor, 1, 1, tracks.size());

	vector<Arr2D> groupResultM, groupResultM_cor;
	AllocGroupMatrices(groupResultM, tracks.size());
	AllocGroupMatrices(groupResultM_cor, tracks.size());

	vector<Beam> outBeams;

	All = Arr2D(1, 1, 4, 4);
	All_cor = Arr2D(1, 1, 4, 4);
	gNorm = gammaRange.step/gammaRange.norm;

	OutputStartTime(timer);

	double beta, gamma;

	for (int i = 0; i <= betaRange.number; ++i)
	{
		m_incomingEnergy = 0;
		OutputProgress(betaRange.number, count, timer);
		++count;

		beta = betaRange.min + betaRange.step*i;

		for (int j = 0; j <= gammaRange.number; ++j)
		{
			gamma = gammaRange.min + gammaRange.step*j;
			m_tracing->SplitBeamByParticle(beta, gamma, outBeams);

			m_incomingEnergy += m_tracing->GetIncomingEnergy();

			HandleBeamsBackScatterPO(outBeams, wave, tracks);
			outBeams.clear();
			OutputState(i, j);

			AddResultToMatrices(groupResultM, gNorm);
			AddResultToMatrix(All, J, gNorm);
			CleanJ(J);

			AddResultToMatrices(groupResultM_cor, gNorm);
			AddResultToMatrix(All_cor, J_cor, gNorm);
			CleanJ(J_cor);
#ifdef _DEBUG // DEB
			cout << "\r     \rj: " << j;
#endif
		}

		m_incomingEnergy *= gNorm;

		double degBeta = RadToDeg(beta);
		OutputToAllFile(diffFile, otherFile, degBeta, allFile, All, Other);
		OutputToAllFile(diffFile_cor, otherFile_cor, degBeta, allFile_cor,
						All_cor, Other_cor);

		OutputToGroupFiles(degBeta, groupFiles, groupResultM, tracks.size());
		OutputToGroupFiles(degBeta, groupFiles_cor, groupResultM_cor, tracks.size());
	}

	allFile.close();
	allFile_cor.close();

	if (isCalcOther)
	{
		otherFile.close();
		diffFile.close();

		otherFile_cor.close();
		diffFile_cor.close();
	}

	for (size_t group = 0; group < groupResultM.size(); ++group)
	{
		ofstream &file = *(groupFiles[group]);
		file.close();

		ofstream &file_cor = *(groupFiles_cor[group]);
		file_cor.close();
	}

	long long orNumber = betaRange.number * gammaRange.number;
	OutputStatisticsPO(timer, orNumber);
}

//REF: объединить с предыдущим
void Tracer::TraceIntervalPO2(int betaNumber, int gammaNumber, const Cone &bsCone,
							  const Tracks &tracks, double wave)
{
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
				HandleBeamsPO2(outBeams, bsCone, wave, groupID);

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

void Tracer::AddResultToMatrices(std::vector<Arr2D> &M, double norm)
{
	for (size_t q = 0; q < J.size(); ++q)
	{
		matrix m = Mueller(J[q](0, 0));
		m *= norm;
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
								 Contribution &contr)
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
	ofstream out("out.dat", ios::out);
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

#ifdef _DEBUG // DEB
//	beta = (155 + 0.5)*betaNorm;
//	gamma = (640 + 0.5)*gammaNorm;
//	m_tracing->SplitBeamByParticle(beta, gamma, outBeams);
#endif

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
			OutputState(i, j);
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
	Arr2D M(bsCone.phiCount+1, bsCone.thetaCount+1, 4, 4);
	ofstream outFile(m_resultDirName, ios::out);
	vector<Beam> outBeams;

	double b = DegToRad(beta);
	double g = DegToRad(gamma);
	m_tracing->SplitBeamByParticle(b, g, outBeams);

	CleanJ(tracks.size(), bsCone);
	HandleBeamsPO(outBeams, bsCone, wave, tracks);

	outBeams.clear();
	AddResultToMatrix(M, bsCone);
	WriteConusMatrices(outFile, M, bsCone);
}

void Tracer::setIsCalcOther(bool value)
{
	isCalcOther = value;
}

void Tracer::SetJnRot(Beam &beam, const Point3f &T,
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
						   double wavelength, const Tracks &tracks)
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
		T = T/Length(T); // базис выходящего пучка

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

				matrixC Jn_rot(2, 2);
				SetJnRot(beam, T, vf, vr, Jn_rot);

				complex fn(0, 0);
				fn = beam.DiffractionIncline(vr, wavelength);

				double dp = DotProductD(vr, Point3d(center));
				complex tmp = exp_im(M_2PI*(lng_proj0-dp)/wavelength);
				matrixC fn_jn = beam.J * tmp;

				matrixC c = fn*Jn_rot*fn_jn;
				J[groupID].insert(i, j, c);
			}
		}
	}
}

void Tracer::HandleBeamsPO2(vector<Beam> &outBeams, const Cone &bsCone,
							double wavelength, int groupID)
{
	for (unsigned int i = 0; i < outBeams.size(); ++i)
	{
		Beam &beam = outBeams.at(i);

		beam.RotateSpherical(-m_incidentDir, m_polarizationBasis);

		Point3f center = beam.Center();
		Point3d center_d = Point3d(center);
		double lng_proj0 = beam.opticalPath + DotProduct(center, beam.direction);

		Point3f T = CrossProduct(beam.e, beam.direction);
		T = T/Length(T); // базис выходящего пучка

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
				SetJnRot(beam, T, vf, vr, Jn_rot);

				complex fn(0, 0);
				fn = beam.DiffractionIncline(vr, wavelength);

				double dp = DotProductD(vr, center_d);
				complex tmp = exp_im(M_2PI*(lng_proj0-dp)/wavelength);
				matrixC fn_jn = beam.J * tmp;

				matrixC c = fn*Jn_rot*fn_jn;
				J[groupID].insert(i, j, c);
			}
		}
	}
}

void Tracer::HandleBeamsBackScatterPO(std::vector<Beam> &outBeams,
									  double wavelength, const Tracks &tracks)
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

		Point3f center = beam.Center();
		double lng_proj0 = beam.opticalPath + DotProduct(center, beam.direction);

		Point3f T = CrossProduct(beam.e, beam.direction);
		T = T/Length(T); // базис выходящего пучка

		matrixC Jn_rot(2, 2);
		SetJnRot(beam, T, vf, vr, Jn_rot);

		complex fn(0, 0);
		fn = beam.DiffractionIncline(vr, wavelength);

		double dp = DotProductD(vr, Point3d(center));
		complex tmp = exp_im(M_2PI*(lng_proj0-dp)/wavelength);
		matrixC fn_jn = beam.J * tmp;

		matrixC c = fn*Jn_rot*fn_jn;

		// correction
		matrixC c_cor = c;
		c_cor[0][1] -= c_cor[1][0];
		c_cor[0][1] /= 2;
		c_cor[1][0] = -c_cor[0][1];

		if (groupID < 0 && isCalcOther)
		{
			matrix m = Mueller(c);
			m *= gNorm;

			Other.insert(0, 0, m);
			All.insert(0, 0, m);

			matrix m_cor = Mueller(c_cor);
			m_cor *= gNorm;

			Other_cor.insert(0, 0, m_cor);
			All_cor.insert(0, 0, m_cor);
		}
		else
		{
			J[groupID].insert(0, 0, c);
			J_cor[groupID].insert(0, 0, c_cor);
		}
	}
}
