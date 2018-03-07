#include "Tracer.h"

#include <ostream>
#include <iostream>
#include <assert.h>
#include "global.h"
#include "macro.h"
//#ifdef _TRACK_ALLOW
//std::ofstream trackMapFile("tracks_deb.dat", std::ios::out);
//#endif

#include "ScatteringConvex.h"
#include "ScatteringNonConvex.h"
#include "ScatteringFiles.h"

using namespace std;

Tracer::Tracer(Particle *particle, int reflNum, const string &resultFileName)
	: m_resultDirName(resultFileName)
{
	SetIncidentLight(particle);

	if (particle->IsConcave())
	{
		m_scattering = new ScatteringNonConvex(particle, &m_incidentLight, true, reflNum);
	}
	else
	{
		m_scattering = new ScatteringConvex(particle, &m_incidentLight, true, reflNum);
	}

	m_particle = m_scattering->m_particle;
	m_symmetry = m_particle->GetSymmetry();
}

Tracer::~Tracer()
{
}

void Tracer::SetIncidentLight(Particle *particle)
{
	m_incidentLight.direction = Point3f(0, 0, -1);
	m_incidentLight.polarizationBasis = Point3f(0, 1, 0);

	Point3f point = m_incidentLight.direction * particle->GetRotationRadius();
	m_incidentLight.direction.d_param = DotProduct(point, m_incidentLight.direction);
}

void Tracer::OutputOrientationToLog(int i, int j, ostream &logfile)
{
	logfile << "i: " << i << ", j: " << j << endl;
	logfile.flush();
}

void Tracer::OutputProgress(int betaNumber, long long count, CalcTimer &timer)
{
	EraseConsoleLine(50);
	cout << (count*100)/(betaNumber+1) << '%'
		 << '\t' << timer.Elapsed();
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
	normIndex = gammaStep/m_symmetry.gamma;

	timer.Start();

	for (int i = 0; i <= betaNumber; ++i)
	{
		beta = i*betaStep;

		for (int j = -halfGammaCount; j <= halfGammaCount; ++j)
		{
			gamma = j*gammaStep;

			m_scattering->ScatterLight(beta, gamma, outBeams);

			CleanJ(tracks.size(), bsCone);
			HandleBeamsPO(outBeams, bsCone, tracks);

			outBeams.clear();
			AddResultToMatrix(M, bsCone, normIndex);
		}

		WriteConusMatrices(outFile, M, bsCone);

		OutputProgress(betaNumber, count, timer);
		++count;
	}

	outFile.close();
}

void Tracer::OutputStatisticsPO(CalcTimer &timer, long long orNumber, const string &path)
{
	string startTime = ctime(&m_startTime);
	string totalTime = timer.Elapsed();
	time_t end = timer.Stop();
	string endTime = ctime(&end);

	m_summary += "\nStart of calculation = " + startTime
			+ "End of calculation   = " + endTime
			+ "\nTotal time of calculation = " + totalTime
			+ "\nTotal number of body orientation = " + to_string(orNumber);

	if (isNanOccured)
	{
		m_summary += "\n\nWARNING! NAN values occured. See 'log.txt'";
	}

	ofstream out(path + "\\out.dat", ios::out);

	out << m_summary;
	out.close();

	cout << m_summary;
}

void Tracer::SetIsOutputGroups(bool value)
{
	isOutputGroups = value;
}

//REF: объединить с предыдущим
void Tracer::TraceRandomPO2(int betaNumber, int gammaNumber, const Cone &bsCone,
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
	normIndex = gammaNorm/m_symmetry.gamma;

	timer.Start();

	for (int i = 0; i <= betaNumber; ++i)
	{
		beta = i*betaNorm;

		for (int j = -halfGammaCount; j <= halfGammaCount; ++j)
		{
			gamma = j*gammaNorm;

			for (size_t groupID = 0; groupID < tracks.size(); ++groupID)
			{
				m_scattering->ScatterLight(beta, gamma, tracks[groupID].tracks, outBeams);

				CleanJ(tracks.size(), bsCone);
				HandleBeamsPO2(outBeams, bsCone, groupID);

				outBeams.clear();
				AddResultToMatrix(M, bsCone, normIndex);
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

void Tracer::OutputStartTime(CalcTimer &timer)
{
	m_startTime = timer.Start();
	cout << "Started at " << ctime(&m_startTime) << endl;
}

void Tracer::TraceFixedPO(const double &beta, const double &gamma,
							 const Cone &bsCone, const Tracks &tracks, double wave)
{
	m_wavelength = wave;
	Arr2D M(bsCone.phiCount+1, bsCone.thetaCount+1, 4, 4);
	ofstream outFile(m_resultDirName, ios::out);
	vector<Beam> outBeams;

	double b = DegToRad(beta);
	double g = DegToRad(gamma);
	m_scattering->ScatterLight(b, g, outBeams);

	CleanJ(tracks.size(), bsCone);
	HandleBeamsPO(outBeams, bsCone, tracks);

	outBeams.clear();
	AddResultToMatrix(M, bsCone);
	WriteConusMatrices(outFile, M, bsCone);
}

void Tracer::SetIsCalcOther(bool value)
{
	isCalcOther = value;
}

void Tracer::CalcJnRot(const Beam &beam, const Point3f &T, const Point3d &vf,
					   const Point3d &vr, matrixC &Jn_rot)
{
	Point3f normal = beam.Normal();

	Point3d vt = CrossProductD(vf, vr);
	vt = vt/LengthD(vt);

	Point3f NT = CrossProduct(normal, T);
	Point3f NE = CrossProduct(normal, beam.polarizationBasis);

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
		int groupID = tracks.FindGroup(beam.trackId);

		if (groupID < 0)
		{
			continue;
		}

		beam.RotateSpherical(-m_incidentLight.direction,
							 m_incidentLight.polarizationBasis);

		Point3f center = beam.Center();
		double lng_proj0 = beam.opticalPath + DotProduct(center, beam.direction);

		Point3f T = CrossProduct(beam.polarizationBasis, beam.direction);
		T = T/Length(T); // basis of beam

		for (int i = 0; i <= bsCone.phiCount; ++i)
		{
			for (int j = 0; j <= bsCone.thetaCount; ++j)
			{	//
				double f = i*bsCone.dPhi;
				double t = j*bsCone.dTheta;

				double sinT = sin(t), sinF = sin(f), cosF = cos(f);

				Point3d vr(sinT*cosF, sinT*sinF, cos(t));
				Point3d vf = (j == 0) ? -m_incidentLight.polarizationBasis
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

		beam.RotateSpherical(-m_incidentLight.direction,
							 m_incidentLight.polarizationBasis);

		Point3f center = beam.Center();
		Point3d center_d = Point3d(center);
		double lng_proj0 = beam.opticalPath + DotProduct(center, beam.direction);

		Point3f T = CrossProduct(beam.polarizationBasis, beam.direction);
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
				Point3d vf = (j == 0) ? -m_incidentLight.polarizationBasis
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
