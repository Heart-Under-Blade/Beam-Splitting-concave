#include "TracerGO.h"
#include <iostream>

using namespace std;

TracerGO::TracerGO(Particle *particle, int maxActNo, const std::string &resultFileName)
	: Tracer(particle, maxActNo, resultFileName)
{
}

void TracerGO::TraceRandom(const AngleRange &betaRange, const AngleRange &gammaRange)
{
#ifdef _CHECK_ENERGY_BALANCE
	m_incomingEnergy = 0;
	m_outcomingEnergy = 0;
#endif
	vector<Beam> scatteredBeams;
	Orientation angle;

	CalcTimer timer;
	OutputStartTime(timer);

	for (int i = 0; i < betaRange.number; ++i)
	{
		angle.beta = (i + 0.5)*betaRange.step;

		for (int j = 0; j < gammaRange.number; ++j)
		{
			angle.gamma = (j + 0.5)*gammaRange.step;
#ifdef _DEBUG // DEB
//			angle.beta = Angle::DegToRad(179.34);
//			angle.gamma = Angle::DegToRad(37);
#endif
			m_particle->Rotate(angle);
			m_scattering->ScatterLight(scatteredBeams);
			m_handler->HandleBeams(scatteredBeams);
			scatteredBeams.clear();

#ifdef _CHECK_ENERGY_BALANCE
			m_incomingEnergy += m_scattering->GetIncidentEnergy()*sin(angle.beta);
#endif
//			m_handler->WriteLog(to_string(i) + ", " + to_string(j) + " ");
//			OutputOrientationToLog(i, j, logfile);
		}

		OutputProgress(betaRange.number, i, timer);
	}

	long long orNum = gammaRange.number * betaRange.number;
	double norm = CalcNorm(orNum);
	m_handler->SetNormIndex(norm);

	m_outcomingEnergy = ((HandlerGO*)m_handler)->ComputeTotalScatteringEnergy();
	m_handler->WriteMatricesToFile(m_resultDirName);
	OutputSummary(orNum, m_outcomingEnergy, norm, timer);
}

void TracerGO::TraceFixed(const double &beta, const double &gamma)
{
	Orientation angle;
	angle.beta = Orientation::DegToRad(beta);
	angle.gamma = Orientation::DegToRad(gamma);

	vector<Beam> outBeams;
	m_particle->Rotate(angle);
	m_scattering->ScatterLight(outBeams);
//	m_particle->Output();
	m_handler->HandleBeams(outBeams);
	outBeams.clear();

//	double D_tot = CalcTotalScatteringEnergy();

	m_handler->WriteMatricesToFile(m_resultDirName);
//	WriteStatisticsToFileGO(1, D_tot, 1, timer); // TODO: раскомментить
}

double TracerGO::CalcNorm(long long orNum)
{
	const double &symBeta = m_particle->GetSymmetry().beta;
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

	m_summary += "\nStart of calculation = " + startTime
			+ "End of calculation   = " + endTime
			+ "\nTotal time of calculation = " + totalTime
			+ "\nTotal number of body orientation = " + to_string(orNumber)
			+ "\nTotal scattering energy = " + to_string(D_tot);

#ifdef _CHECK_ENERGY_BALANCE
	double normEnergy = m_incomingEnergy * NRM;
	double passedEnergy = (m_outcomingEnergy/normEnergy)*100;

	m_summary += "\nTotal incoming energy = " + to_string(normEnergy)
			+ "\nTotal outcoming energy = " + to_string(m_outcomingEnergy)
			+ "\nEnergy passed = " + to_string(passedEnergy) + '%';
#endif

	// out << "\nAveraged cross section = " << incomingEnergy*NRM;
	ofstream out(m_resultDirName+"_out.dat", ios::out);
	out << m_summary;
	out.close();

	cout << m_summary;
}
