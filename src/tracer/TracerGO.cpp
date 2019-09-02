#include "TracerGO.h"
<<<<<<< HEAD
#include "handler/HandlerGO.h"
=======

#include "HandlerGO.h"

>>>>>>> origin/refactor
#include <iostream>

using namespace std;

<<<<<<< HEAD
TracerGO::TracerGO(Particle *particle, int reflNum, const std::string &resultFileName)
	: Tracer(particle, reflNum, resultFileName)
{
}

void TracerGO::TraceRandom(const AngleRange &betaRange, const AngleRange &gammaRange)
=======
TracerGO::TracerGO(Particle *particle, Scattering *scattering,
				   const std::string &resultFileName)
	: Tracer(particle, scattering, resultFileName)
{
}

void TracerGO::TraceRandom(const OrientationRange &range)
>>>>>>> origin/refactor
{
#ifdef _CHECK_ENERGY_BALANCE
	m_incomingEnergy = 0;
	m_outcomingEnergy = 0;
#endif
<<<<<<< HEAD

	vector<Beam> outBeams;
	double beta, gamma;
=======
	long long orNum = range.nZenith * range.nAzimuth;

	vector<Beam> beams;
	Orientation orientation;
>>>>>>> origin/refactor

	CalcTimer timer;
	OutputStartTime(timer);

<<<<<<< HEAD
	for (int i = 0; i < betaRange.number; ++i)
	{
		beta = (i + 0.5)*betaRange.step;

		for (int j = 0; j < gammaRange.number; ++j)
		{
			gamma = (j + 0.5)*gammaRange.step;

			m_particle->Rotate(beta, gamma, 0);
			m_scattering->ScatterLight(outBeams);
			m_handler->HandleBeams(outBeams);
			outBeams.clear();
=======
	long long count = 0;

	for (int i = 0; i < range.nZenith; ++i)
	{
		orientation.zenith = (i + 0.5)*range.step.zenith;

		for (int j = 0; j < range.nAzimuth; ++j)
		{
			orientation.azimuth = (j + 0.5)*range.step.azimuth;
#ifdef _DEBUG // DEB
//			angle.beta = Angle::DegToRad(179.34);
//			angle.gamma = Angle::DegToRad(37);
#endif
			m_particle->Rotate(orientation);
			m_scattering->ScatterLight(beams);
			m_handler->HandleBeams(beams);
			beams.clear();
>>>>>>> origin/refactor

#ifdef _CHECK_ENERGY_BALANCE
			m_incomingEnergy += m_scattering->GetIncedentEnergy()*sin(beta);
#endif
//			m_handler->WriteLog(to_string(i) + ", " + to_string(j) + " ");
//			OutputOrientationToLog(i, j, logfile);
			OutputProgress(orNum, count++, timer);
		}
<<<<<<< HEAD

		OutputProgress(betaRange.number, i, timer);
	}

	long long orNum = gammaRange.number * betaRange.number;
=======
	}

>>>>>>> origin/refactor
	double norm = CalcNorm(orNum);
	m_handler->SetNormIndex(norm);

	m_outcomingEnergy = ((HandlerGO*)m_handler)->ComputeTotalScatteringEnergy();
	m_handler->WriteMatricesToFile(m_resultDirName);
	OutputSummary(orNum, m_outcomingEnergy, norm, timer);
}

void TracerGO::TraceFixed(const double &beta, const double &gamma)
{
	double b = DegToRad(beta);
	double g = DegToRad(gamma);

	vector<Beam> outBeams;
	m_particle->Rotate(b, g, 0);
	m_scattering->ScatterLight(outBeams);
	m_handler->HandleBeams(outBeams);
	outBeams.clear();

//	double D_tot = CalcTotalScatteringEnergy();

	m_handler->WriteMatricesToFile(m_resultDirName);
//	WriteStatisticsToFileGO(1, D_tot, 1, timer); // TODO: раскомментить
}

double TracerGO::CalcNorm(long long orNum)
{
	double &symBeta = m_symmetry.beta;
	double tmp = (/*isRandom*/true) ? symBeta : 1.0;
	double dBeta = -(cos(symBeta) - cos(0));
	return tmp/(orNum*dBeta);
}

void TracerGO::OutputSummary(int orNumber, double D_tot, double NRM, CalcTimer &timer)
{
	string startTime = ctime(&m_startTime);
	string totalTime = timer.Elapsed();
	time_t end = timer.Stop();
	string endTime = ctime(&end);

	m_log += "\nStart of calculation = " + startTime
			+ "End of calculation   = " + endTime
			+ "\nTotal time of calculation = " + totalTime
			+ "\nTotal number of body orientation = " + to_string(orNumber)
			+ "\nTotal scattering energy = " + to_string(D_tot);

#ifdef _CHECK_ENERGY_BALANCE
	double normEnergy = m_incomingEnergy * NRM;
	double passedEnergy = (m_outcomingEnergy/normEnergy)*100;

	m_log += "\nTotal incoming energy = " + to_string(normEnergy)
			+ "\nTotal outcoming energy = " + to_string(m_outcomingEnergy)
			+ "\nEnergy passed = " + to_string(passedEnergy) + '%';
#endif

	// out << "\nAveraged cross section = " << incomingEnergy*NRM;
	ofstream out(m_resultDirName+"_out.dat", ios::out);
	out << m_log;
	out.close();

	cout << m_log;
}
