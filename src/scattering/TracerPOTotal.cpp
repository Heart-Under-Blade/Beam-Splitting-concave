#include "TracerPOTotal.h"

#include <iostream>

using namespace std;

TracerPOTotal::TracerPOTotal(Particle *particle, int nActs,
							 const string &resultFileName)
	: TracerPO(particle, nActs, resultFileName)
{
}

void TracerPOTotal::TraceRandom(const AngleRange &betaRange,
								const AngleRange &gammaRange)
{
	CalcTimer timer;
	long long count = 0;
	long long nOrientations = betaRange.number * gammaRange.number;

	ofstream outFile(m_resultDirName + ".dat", ios::out);

	if (!outFile.is_open())
	{
		std::cerr << "Error! File \"" << m_resultDirName
				  << "\" was not opened. " << __FUNCTION__;
		throw std::exception();
	}

	vector<Beam> outBeams;
	double beta, gamma;

	double normIndex = 2 * gammaRange.number;
	m_handler->SetNormIndex(normIndex);

	++nOrientations;
	timer.Start();

	for (int i = 0; i <= betaRange.number; ++i)
	{
		beta = i*betaRange.step;

		double sinZenith = (i == 0 || i == betaRange.number)
				? (1.0-cos(0.5*betaRange.step))/normIndex
				: (cos((i-0.5)*betaRange.step) -
				   cos((i+0.5)*betaRange.step))/normIndex;

		m_handler->SetSinZenith(sinZenith);

		for (int j = 0; j < gammaRange.number; ++j)
		{
			gamma = j*gammaRange.step;

			m_particle->Rotate(M_PI-beta, M_PI+gamma, 0);
			m_scattering->FormShadowBeam(outBeams);
			m_scattering->ScatterLight(M_PI-beta, M_PI+gamma, outBeams);

			m_handler->HandleBeams(outBeams);
			outBeams.clear();

			OutputProgress(nOrientations, count, beta, gamma, timer);
			++count;
		}
	}

	std::string dir = CreateFolder(m_resultDirName);
	m_resultDirName = dir + m_resultDirName + '\\' + m_resultDirName;
	m_handler->WriteMatricesToFile(m_resultDirName);
	OutputStatisticsPO(timer, nOrientations, dir);
	outFile.close();
}
