#include "TracerGO.h"
#include <iostream>

using namespace std;

TracerGO::TracerGO(Particle *particle, int reflNum, const std::string &resultFileName)
	: Tracer(particle, reflNum, resultFileName)
{
}

void TracerGO::TraceRandom(const AngleRange &betaRange, const AngleRange &gammaRange,
						   bool isCalcTracks)
{
#ifdef _CHECK_ENERGY_BALANCE
	m_incomingEnergy = 0;
	m_outcomingEnergy = 0;
#endif

	vector<Beam> outBeams;
	double beta, gamma;

	SetHandler(isCalcTracks);

	CalcTimer timer;
	OutputStartTime(timer);

	for (int i = 0; i < betaRange.number; ++i)
	{
		beta = (i + 0.5)*betaRange.step;

		for (int j = 0; j < gammaRange.number; ++j)
		{
			gamma = (j + 0.5)*gammaRange.step;
			m_tracing->SplitBeamByParticle(beta, gamma, outBeams);
			m_handler->HandleBeams(outBeams, beta);

#ifdef _CHECK_ENERGY_BALANCE
			m_incomingEnergy += m_tracing->GetIncomingEnergy()*sin(beta);
#endif
//			OutputOrientationToLog(i, j, logfile);
		}

		OutputProgress(betaRange.number, i, timer);
	}

	long long orNum = gammaRange.number * betaRange.number;
	double norm = CalcNorm(orNum);

	m_outcomingEnergy = m_handler->CalcTotalScatteringEnergy(norm);
	m_handler->WriteMatricesToFile(norm, m_resultDirName);
	OutputSummary(orNum, m_outcomingEnergy, norm, timer);
}

void TracerGO::TraceFixed(const double &beta, const double &gamma, bool isCalcTracks)
{
	SetHandler(isCalcTracks);

	double b = DegToRad(beta);
	double g = DegToRad(gamma);

	vector<Beam> outBeams;
	m_tracing->SplitBeamByParticle(b, g, outBeams);
	m_handler->HandleBeams(outBeams, beta);

//	double D_tot = CalcTotalScatteringEnergy();

	m_handler->WriteMatricesToFile(1, m_resultDirName);
//	WriteStatisticsToFileGO(1, D_tot, 1, timer); // TODO: раскомментить
}

void TracerGO::SetTracks(Tracks *tracks)
{
	m_tracks = tracks;
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

	//	out << "\nAveraged cross section = " << incomingEnergy*NRM;
	ofstream out(m_resultDirName+"_out.dat", ios::out);
	out << m_summary;
	out.close();

	cout << m_summary;
}

void TracerGO::SetHandler(bool isCalcTracks)
{
	if (isCalcTracks)
	{
		m_handler = new HandlerGroupGO(m_tracing->m_particle, &m_incidentLight);
		m_handler->SetTracks(m_tracks);
	}
	else
	{
		m_handler = new HandlerTotalGO(m_tracing->m_particle, &m_incidentLight);
	}
}
