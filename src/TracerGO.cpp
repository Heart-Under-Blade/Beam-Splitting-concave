#include "TracerGO.h"
#include <iostream>

using namespace std;

TracerGO::TracerGO(Particle *particle, int reflNum, const std::string &resultFileName)
	: Tracer(particle, reflNum, resultFileName)
{
}

void TracerGO::TraceRandom(const AngleRange &betaRange, const AngleRange &gammaRange)
{
#ifdef _CHECK_ENERGY_BALANCE
	m_incomingEnergy = 0;
	m_outcomingEnergy = 0;
#endif
	vector<Beam> outBeams;
	double beta, gamma;

	CalcTimer timer;
	OutputStartTime(timer);

	for (int i = 0; i < betaRange.number; ++i)
	{
		beta = (i + 0.5)*betaRange.step;

		for (int j = 0; j < gammaRange.number; ++j)
		{
			gamma = (j + 0.5)*gammaRange.step;
#ifdef _DEBUG // DEB
//			beta = 0.47123889803846897; gamma = 0.52359877559829882;
#endif
			m_scattering->ScatterLight(beta, gamma, outBeams);
//			m_particle->Output();
			m_handler->HandleBeams(outBeams);
			outBeams.clear();

#ifdef _CHECK_ENERGY_BALANCE
			m_incomingEnergy += m_scattering->GetIncedentEnergy()*sin(beta);
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

	std::string dir = CreateFolder(m_resultDirName);
	m_resultDirName = dir + m_resultDirName + '\\' + m_resultDirName;

	m_handler->WriteMatricesToFile(m_resultDirName);
	OutputSummary(orNum, m_outcomingEnergy, norm, timer);
}

void TracerGO::TraceFixed(const double &beta, const double &gamma)
{
	double b = DegToRad(beta);
	double g = DegToRad(gamma);

	vector<Beam> outBeams;
	m_scattering->ScatterLight(b, g, outBeams);
	m_handler->HandleBeams(outBeams);
	outBeams.clear();

//	double D_tot = CalcTotalScatteringEnergy();

	m_handler->WriteMatricesToFile(m_resultDirName);
//	WriteStatisticsToFileGO(1, D_tot, 1, timer); // TODO: раскомментить
}

double TracerGO::CalcNorm(long long orNum)
{
//	return 1.0/(double)orNum;

	double &symBeta = m_symmetry.beta;

//	if (symBeta < M_PI - FLT_EPSILON) //
	{
		double tmp = symBeta;
		double dBeta = -(cos(symBeta) - 1);
//		dBeta = (dBeta < 0) ? -dBeta : dBeta;
		return tmp/(orNum*dBeta);
	}
//	else // otherwise the result becomes 'inf'
	{
//		return 1;
	}
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
	const double normEnergy = m_incomingEnergy * NRM;
	const double passedEnergy = (m_outcomingEnergy/normEnergy)*100;
	const double parArea = m_particle->Area()/4.0;

	m_summary += "\nTotal incoming energy = " + to_string(normEnergy)
			+ " must be equal to " + to_string(parArea) + " (S/4)"
			+ "\nTotal outcoming energy = " + to_string(m_outcomingEnergy)
			+ "\nEnergy passed = " + to_string(passedEnergy) + '%'
			+ "\nNorm index = " + to_string(NRM);
#endif

	// out << "\nAveraged cross section = " << incomingEnergy*NRM;
	ofstream out(m_resultDirName+"_out.dat", ios::out);
	out << m_summary;
	out.close();

	cout << m_summary;
}
