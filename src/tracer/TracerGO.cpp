#include "TracerGO.h"
#include "HandlerGO.h"
#include <iostream>

using namespace std;

TracerGO::TracerGO(Particle *particle, Scattering *scattering,
				   const std::string &resultFileName)
	: LightTracer(particle, scattering, resultFileName)
{
}

void TracerGO::TraceRandom(const AngleRange &zenithRange,
						   const AngleRange &azimuthRange)
{
#ifdef _CHECK_ENERGY_BALANCE
	m_incomingEnergy = 0;
	m_outcomingEnergy = 0;
#endif
	vector<Beam> beams;
	Orientation orientation;

	CalcTimer timer;
	OutputStartTime(timer);

	for (int i = 0; i < zenithRange.number; ++i)
	{
		orientation.zenith = (i + 0.5)*zenithRange.step;

		for (int j = 0; j < azimuthRange.number; ++j)
		{
			orientation.azimuth = (j + 0.5)*azimuthRange.step;
#ifdef _DEBUG // DEB
//			angle.beta = Angle::DegToRad(179.34);
//			angle.gamma = Angle::DegToRad(37);
#endif
			m_particle->Rotate(orientation);
			m_scattering->ScatterLight(beams);
			m_handler->HandleBeams(beams);
			beams.clear();

#ifdef _CHECK_ENERGY_BALANCE
			m_incomingEnergy += m_scattering->GetIncidentEnergy()*sin(orientation.zenith);
#endif
//			m_handler->WriteLog(to_string(i) + ", " + to_string(j) + " ");
//			OutputOrientationToLog(i, j, logfile);
		}

		OutputProgress(zenithRange.number, i, timer);
	}

	long long orNum = azimuthRange.number * zenithRange.number;
	double norm = CalcNorm(orNum);
	m_handler->SetNormIndex(norm);

	m_outcomingEnergy = ((HandlerGO*)m_handler)->ComputeTotalScatteringEnergy();
	m_handler->WriteMatricesToFile(m_resultDirName);
	OutputSummary(orNum, m_outcomingEnergy, norm, timer);
}

double TracerGO::CalcNorm(long long orNum)
{
	const double &symZenith = m_particle->GetSymmetry().zenith;
	double tmp = (/*isRandom*/true) ? symZenith : 1.0;
	double dBeta = -(cos(symZenith) - cos(0));
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
