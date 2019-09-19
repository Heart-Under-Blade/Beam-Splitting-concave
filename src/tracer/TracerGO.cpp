#include "TracerGO.h"

#include <iostream>

#include "HandlerGO.h"

using namespace std;

TracerGO::TracerGO(Particle *particle, Scattering *scattering,
				   const std::string &resultFileName)
	: Tracer(particle, scattering, resultFileName)
{
}

void TracerGO::TraceRandom(const OrientationRange &range)
{
#ifdef _CHECK_ENERGY_BALANCE
	m_incomingEnergy = 0;
	m_outcomingEnergy = 0;
#endif

	vector<Beam> outBeams;
	double beta, gamma;

	long long orNum = range.nZenith * range.nAzimuth;

	vector<Beam> beams;
	Orientation orientation;

	CalcTimer timer;
	OutputStartTime(timer);

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
			m_scattering->Scatter(&beams);
			m_handler->HandleBeams(beams);
			beams.clear();

#ifdef _CHECK_ENERGY_BALANCE
			m_incomingEnergy += m_scattering->GetIncidentEnergy()*sin(beta);
#endif
//			m_handler->WriteLog(to_string(i) + ", " + to_string(j) + " ");
//			OutputOrientationToLog(i, j, logfile);
			OutputProgress(orNum, count++, timer);
		}

		OutputProgress(range.nZenith, i, timer);
	}

	double norm = CalcNorm(orNum);
	m_handler->SetNormIndex(norm);

	m_outcomingEnergy = ((HandlerGO*)m_handler)->ComputeTotalScatteringEnergy();
	m_handler->WriteMatricesToFile(m_resultDirName);
	OutputSummary(orNum, m_outcomingEnergy, norm, timer);
}

void TracerGO::TraceFixed(const Orientation &orientation)
{
	double b = Orientation::DegToRad(orientation.zenith);
	double g = Orientation::DegToRad(orientation.azimuth);

	vector<Beam> outBeams;
	m_particle->Rotate(orientation);
	m_scattering->Scatter(&outBeams);
	m_handler->HandleBeams(outBeams);
	outBeams.clear();

//	double D_tot = CalcTotalScatteringEnergy();

	m_handler->WriteMatricesToFile(m_resultDirName);
//	WriteStatisticsToFileGO(1, D_tot, 1, timer); // TODO: раскомментить
}

double TracerGO::CalcNorm(long long orNum)
{
	double &symBeta = m_symmetry.zenith;
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
