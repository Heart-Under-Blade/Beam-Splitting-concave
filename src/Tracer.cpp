#include "Tracer.h"

#include <iostream>

using namespace std;

Tracer::Tracer(Tracing *tracing, const string resultFileName)
	: m_tracing(tracing),
	  m_incidentDir(0, 0, -1),
	  m_polarizationBasis(0, 1, 0),
	  m_resultFileName(resultFileName),
	  m_gammaNorm(tracing->m_particle->GetSymmetryGamma())
{
}

void Tracer::WriteResultToSeparateFilesGO(double NRM, int thetaNum, int EDF,
										  const Tracks &tracks)
{
	for (int i = 0; i < m_sepatateMatrices.size(); ++i)
	{
		if (tracks[i].size != 0)
		{
			string subname;

			if (tracks[i].size == 1)
			{
				for (int index : tracks[i].tracks[0])
				{
					subname += '_' + to_string(index);
				}

				subname += "__" + to_string(i);
			}
			else
			{
				subname += "_group_" + to_string(i);
			}

			ExtractPeaksGO(EDF, NRM, thetaNum, m_sepatateMatrices[i]);
			WriteResultsToFileGO(thetaNum, NRM, m_resultFileName + subname,
								 m_sepatateMatrices[i]);
		}
	}
}

void Tracer::TraceIntervalGO(const AngleRange &betaR, const AngleRange &gammaR,
							 int thetaNum, const Tracks &tracks)
{
	int EDF = 0;
	CalcTimer timer;

#ifdef _CHECK_ENERGY_BALANCE
	m_incomingEnergy = 0;
	m_outcomingEnergy = 0;
#endif

	m_sepatateMatrices.resize(tracks.GetMaxGroupID());

	sizeBin = M_PI/thetaNum;

//	m_totalMtrx.scatMatrix = Arr2D(1, thetaNum+1, 4, 4);
//	m_totalMtrx.scatMatrix.ClearArr();

	vector<Beam> outBeams;
	double beta, gamma;

	m_startTime = timer.Start();
	cout << "Started at " << ctime(&m_startTime) << endl;

	for (int i = 0; i < betaR.count; ++i)
	{
		beta = (i + 0.5)*betaR.norm;

		for (int j = 0; j < gammaR.count; ++j)
		{
			gamma = (j + 0.5)*gammaR.norm;
			m_tracing->SplitBeamByParticle(beta, gamma, outBeams);

#ifdef _CHECK_ENERGY_BALANCE
			m_incomingEnergy += m_tracing->GetIncomingEnergy()*sin(beta);
#endif
			HandleBeamsGO(outBeams, beta, tracks);
			outBeams.clear();
		}

		PrintProgress(betaR.count, i, timer);
	}

	double D_tot = CalcTotalScatteringEnergy(thetaNum);
	long long orNum = gammaR.count * betaR.count;
	double NRM = CalcNorm(orNum);

#ifdef _CHECK_ENERGY_BALANCE
	for (int deg = 0; deg <= thetaNum; ++deg)
	{
		m_outcomingEnergy += m_muller(0, deg, 0, 0)*NRM;
	}
#endif

	WriteResultToSeparateFilesGO(NRM, thetaNum, EDF, tracks);

	ExtractPeaksGO(EDF, NRM, thetaNum, m_totalMtrx);
	WriteResultsToFileGO(thetaNum, NRM, m_resultFileName + "_all", m_totalMtrx);

	WriteStatisticsToFileGO(orNum, D_tot, NRM);
//	WriteStatisticsToConsole(orNum, D_tot, NRM);
}

void Tracer::PrintProgress(int betaNumber, long long count, CalcTimer &timer)
{
	EraseConsoleLine(50);
	if ((count*100)/betaNumber == 99)
		int fff = 0;
	cout << (count*100)/betaNumber << '%'
		 << '\t' << timer.Elapsed();
}

void Tracer::ExtractPeaksGO(int EDF, double NRM, int ThetaNumber,
							Contribution &contr)
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
		contr.scatMatrix(0,ThetaNumber,0,0) += f[0];
		contr.scatMatrix(0,0,0,0) += b[0];
		contr.scatMatrix(0,ThetaNumber,1,1) += f[1];
		contr.scatMatrix(0,0,1,1) += b[1];
		contr.scatMatrix(0,ThetaNumber,2,2) += f[1];
		contr.scatMatrix(0,0,2,2) -= b[1];
		contr.scatMatrix(0,ThetaNumber,3,3) += f[2];
		contr.scatMatrix(0,0,3,3) += b[2];
	}
}

void Tracer::WriteResultsToFileGO(int thetaNum, double NRM, const string &filename,
								  Contribution &contr)
{
	string name = GetFileName(filename);
	ofstream M(name, std::ios::out);

	M << "tetta M11 M12/M11 M21/M11 M22/M11 M33/M11 M34/M11 M43/M11 M44/M11";

	for (int j = thetaNum; j >= 0; j--)
	{
		double sn;

		//Special case in first and last step
		M << '\n' << 180.0/thetaNum*(thetaNum-j) + (j==0 ?-0.25*180.0/thetaNum:0)+(j==(int)thetaNum ?0.25*180.0/thetaNum:0);
		sn = (j==0 || j==(int)thetaNum) ? 1-cos(sizeBin/2.0)
										: (cos((j-0.5)*sizeBin)-cos((j+0.5)*sizeBin));
		matrix bf = contr.scatMatrix(0, j);

		if(bf[0][0] <= DBL_EPSILON)
		{
			M << " 0 0 0 0 0 0 0 0";
		}
		else
		{
			M << ' ' << bf[0][0]*NRM/(2.0*M_PI*sn)
					<< ' ' << bf[0][1]/bf[0][0]
					<< ' ' << bf[1][0]/bf[0][0]
					<< ' ' << bf[1][1]/bf[0][0]
					<< ' ' << bf[2][2]/bf[0][0]
					<< ' ' << bf[2][3]/bf[0][0]
					<< ' ' << bf[3][2]/bf[0][0]
					<< ' ' << bf[3][3]/bf[0][0];
		}
	}

	M.close();
}

void Tracer::WriteStatisticsToFileGO(int orNumber, double D_tot, double NRM)
{
	ofstream out("out.dat", ios::out);

	// Information for log-file
	out << "\nStart of calculation = " << ctime(&m_startTime);
//	out << "\nTotal time of calculation = " << t/CLOCKS_PER_SEC << " seconds";
	out << "\nTotal number of body orientation = " << orNumber;
	out << "\nTotal scattering energy = " << D_tot;

#ifdef _CHECK_ENERGY_BALANCE
	out << "\nTotal incoming energy = " << (m_incomingEnergy *= NRM);
	out << "\nTotal outcoming energy = " << m_outcomingEnergy;
	out << "\nEnergy passed = " << (m_outcomingEnergy/m_incomingEnergy)*100 << '%';
#endif

//	out << "\nAveraged cross section = " << incomingEnergy*NRM;
	out.close();
}

string Tracer::GetFileName(const string &filename)
{
	string fname = string("M_") + filename;
	string name = fname;

	for (int i = 1; ifstream(name += ".dat") != NULL; ++i)
	{
		name = fname + "(" + to_string(i) + ")";
	}

	return name;
}

void Tracer::TraceIntervalPO(const AngleRange &betaR, const AngleRange &gammaR,
							 const Cone &bsCone, const Tracks &tracks, double wave)
{
	CalcTimer timer;
	long long count = 0;

	Arr2D M(bsCone.phiCount+1, bsCone.thetaCount+1, 4, 4);
	ofstream outFile(m_resultFileName, ios::out);

	vector<Beam> outBeams;
	double beta, gamma;
	double bStep = betaR.step;

	int maxGroupID = tracks.GetMaxGroupID();
	int halfGammaCount = gammaR.count/2;
	double gNorm = gammaR.norm*m_gammaNorm;

	timer.Start();

	for (int i = 0; i <= betaR.count; ++i)
	{
		beta = bStep*i;

		for (int j = -halfGammaCount; j <= halfGammaCount; ++j)
		{
			gamma = j*gammaR.norm + M_PI/6;
//beta = DegToRad(32); gamma = DegToRad(30); // DEB
//cout << j << endl;// DEB
			m_tracing->SplitBeamByParticle(beta, gamma, outBeams);

			CleanJ(maxGroupID, bsCone);
			HandleBeamsPO(outBeams, bsCone, wave, tracks);

			outBeams.clear();
			AddResultToSumMatrix(M, maxGroupID, bsCone, gNorm);
		}

		WriteSumMatrix(outFile, M, bsCone);

		PrintProgress(betaR.count, count, timer);
		++count;
	}

	outFile.close();
}

//REF: объединить с предыдущим
void Tracer::TraceIntervalPO2(const AngleRange &betaR, const AngleRange &gammaR,
							  const Cone &bsCone, const Tracks &tracks, double wave)
{
	CalcTimer timer;
	long long count = 0;

	Arr2D M(bsCone.phiCount+1, bsCone.thetaCount+1, 4, 4);
	ofstream outFile(m_resultFileName, ios::out);

	vector<Beam> outBeams;
	double beta, gamma;
	double bStep = betaR.step;

	int maxGroupID = tracks.GetMaxGroupID();
	int halfGammaCount = gammaR.count/2;
	double gNorm = gammaR.norm*m_gammaNorm;

	timer.Start();

	for (int i = 0; i <= betaR.count; ++i)
	{
		beta = bStep*i;
//cout << endl;// DEB
		for (int j = -halfGammaCount; j <= halfGammaCount; ++j)
		{
			gamma = j*gammaR.norm + M_PI/6;

//EraseConsoleLine(50); // DEB
//cout << j;
			for (int groupID = 0; groupID < maxGroupID; ++groupID)
			{
				m_tracing->SplitBeamByParticle(beta, gamma, tracks[groupID].tracks, outBeams);

				CleanJ(maxGroupID, bsCone);
				HandleBeamsPO2(outBeams, bsCone, wave, groupID);

				outBeams.clear();
				AddResultToSumMatrix(M, maxGroupID, bsCone, gNorm);
			}
		}

		WriteSumMatrix(outFile, M, bsCone);

		PrintProgress(betaR.count, count, timer);
		++count;
	}

	outFile.close();
}

void Tracer::WriteSumMatrix(ofstream &outFile, const Arr2D &sum,
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

void Tracer::AddResultToSumMatrix(Arr2D &M_, int maxGroupID, const Cone &bsCone,
								  double norm)
{
	for (int q = 0; q < maxGroupID; ++q)
	{
		for (int t = 0; t <= bsCone.thetaCount; ++t)
		{
			for (int p = 0; p <= bsCone.phiCount; ++p)
			{
//				complex ee = J[q](p, t)[0][0]; // DEB
				matrix Mk = Mueller(J[q](p, t));
				M_.insert(p, t, norm*Mk);
			}
		}
	}
}

void Tracer::CleanJ(int maxGroupID, const Cone &bsCone)
{
	J.clear();
	Arr2DC tmp(bsCone.phiCount+1, bsCone.thetaCount+1, 2, 2);
	tmp.ClearArr();

	for(int q = 0; q < maxGroupID; q++)
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
		const unsigned int zenAng = round(acos(z)/sizeBin);
		contr.scatMatrix.insert(0, zenAng, area*bf);
	}
}

void Tracer::HandleBeamsGO(std::vector<Beam> &outBeams, double beta,
						   const Tracks &tracks)
{
	double sinBeta = sin(beta);

	for (Beam &beam : outBeams)
	{
		int grID = tracks.GetGroupID(beam.id);

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
	double symBeta = m_tracing->m_particle->GetSymmetryBeta();
	double tmp = (/*isRandom*/true) ? symBeta : 1.0;
	double dBeta = -(cos(symBeta) - cos(0));
	return tmp/(orNum*dBeta);
}

double Tracer::CalcTotalScatteringEnergy(int thetaNum)
{
	double D_tot = m_totalMtrx.back[0][0] + m_totalMtrx.forw[0][0];

	for (int i = 0; i <= thetaNum; ++i)
	{
		D_tot += m_totalMtrx.scatMatrix(0, i, 0, 0);
	}

	return D_tot;
}

void Tracer::TraceIntervalGO(const AngleRange &betaR, const AngleRange &gammaR,
							 int thetaNum)
{
	int EDF = 0;
	CalcTimer timer;

#ifdef _CHECK_ENERGY_BALANCE
	m_incomingEnergy = 0;
	m_outcomingEnergy = 0;
#endif

	m_totalMtrx.back.Fill(0);
	m_totalMtrx.forw.Fill(0);

	sizeBin = M_PI/thetaNum;
	m_totalMtrx.scatMatrix = Arr2D(1, thetaNum+1, 4, 4);
	m_totalMtrx.scatMatrix.ClearArr();

	vector<Beam> outBeams;
	double beta, gamma;

	m_startTime = timer.Start();
	cout << "Started at " << ctime(&m_startTime) << endl;

	for (int i = 0; i < betaR.count; ++i)
	{
		beta = (i + 0.5)*betaR.norm;

		for (int j = 0; j < gammaR.count; ++j)
		{
			gamma = (j + 0.5)*gammaR.norm;
			m_tracing->SplitBeamByParticle(beta, gamma, outBeams);

#ifdef _CHECK_ENERGY_BALANCE
			m_incomingEnergy += m_tracing->GetIncomingEnergy()*sin(beta);
#endif
			HandleBeamsGO(outBeams, beta);
			outBeams.clear();
		}

		PrintProgress(betaR.count, i, timer);
	}

	double D_tot = CalcTotalScatteringEnergy(thetaNum);
	long long orNum = gammaR.count * betaR.count;
	double NRM = CalcNorm(orNum);

#ifdef _CHECK_ENERGY_BALANCE
	for (int deg = 0; deg <= thetaNum; ++deg)
	{
		m_outcomingEnergy += m_muller(0, deg, 0, 0)*NRM;
	}
#endif

	ExtractPeaksGO(EDF, NRM, thetaNum, m_totalMtrx);

	WriteResultsToFileGO(thetaNum, NRM, m_resultFileName + "_all", m_totalMtrx);
	WriteStatisticsToFileGO(orNum, D_tot, NRM);
//	WriteStatisticsToConsole(orNum, D_tot, NRM);
}

void Tracer::TraceSingleOrGO(const double &beta, const double &gamma,
							 int thetaNum, const Tracks &tracks)
{
	int EDF = 0;
	vector<Beam> outBeams;

	m_totalMtrx.back.Fill(0);
	m_totalMtrx.forw.Fill(0);

	sizeBin = M_PI/thetaNum;
	m_totalMtrx.scatMatrix = Arr2D(1, thetaNum+1, 4, 4);
	m_totalMtrx.scatMatrix.ClearArr();

	double b = DegToRad(beta);
	double g = DegToRad(gamma);
	m_tracing->SplitBeamByParticle(b, g, outBeams);

	HandleBeamsGO(outBeams, beta, tracks);

	double D_tot = CalcTotalScatteringEnergy(thetaNum);

	ExtractPeaksGO(EDF, 1, thetaNum, m_totalMtrx);

	WriteResultsToFileGO(thetaNum, 1, m_resultFileName, m_totalMtrx);
	WriteStatisticsToFileGO(1, D_tot, 1);
}

void Tracer::TraceSingleOrPO(const double &beta, const double &gamma,
							 const Cone &bsCone, const Tracks &tracks, double wave)
{
	Arr2D M(bsCone.phiCount+1, bsCone.thetaCount+1, 4, 4);
	ofstream outFile(m_resultFileName, ios::out);
	vector<Beam> outBeams;

	double b = DegToRad(beta);
	double g = DegToRad(gamma);
	m_tracing->SplitBeamByParticle(b, g, outBeams);

	int maxGroupID = tracks.GetMaxGroupID();
	CleanJ(maxGroupID, bsCone);
	HandleBeamsPO(outBeams, bsCone, wave, tracks);

	outBeams.clear();
	AddResultToSumMatrix(M, maxGroupID, bsCone);
	WriteSumMatrix(outFile, M, bsCone);
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
//		double ctetta = DotProduct(beam.direction, -m_incidentDir);

//		if (ctetta < 0.17364817766693034885171662676931)
//		{	// отбрасываем пучки, которые далеко от конуса направления назад
//			continue;// DEB
//		}

		int groupID = tracks.GetGroupID(beam.id);

		if (groupID < 0)
		{
			continue;
		}

//		double area = beam.polygon.Area(); // DEB
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
