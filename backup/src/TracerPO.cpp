#include "TracerPO.h"

using namespace std;

TracerPO::TracerPO(Particle *particle, int nActs, const string &resultFileName)
	: Tracer(particle, nActs, resultFileName)
{
}

void TracerPO::TraceRandom(const AngleRange &betaRange, const AngleRange &gammaRange)
{
	CalcTimer timer;
	long long count = 0;

	ofstream outFile(m_resultDirName, ios::out);

	vector<Beam> outBeams;
	double beta, gamma;
	int halfGammaNum = gammaRange.number/2;

	double normIndex = gammaRange.step/gammaRange.norm;
	m_handler->SetNormIndex(normIndex);

	timer.Start();

	for (int i = 0; i <= betaRange.number; ++i)
	{
		beta = i*betaRange.step;

		for (int j = -halfGammaNum; j <= halfGammaNum; ++j)
		{
			gamma = j*gammaRange.step;

			m_scattering->ScatterLight(beta, gamma, outBeams);

			m_handler->HandleBeams(outBeams);
			outBeams.clear();
		}

		m_handler->WriteMatricesToFile(m_resultDirName);

		OutputProgress(betaRange.number, count, timer);
		++count;
	}

	outFile.close();
}

void TracerPO::TraceFixed(const double &beta, const double &gamma)
{
	ofstream outFile(m_resultDirName, ios::out);
	vector<Beam> outBeams;

	double b = DegToRad(beta);
	double g = DegToRad(gamma);
	m_scattering->ScatterLight(b, g, outBeams);

	m_handler->HandleBeams(outBeams);
	outBeams.clear();
	m_handler->WriteMatricesToFile(m_resultDirName);
	outFile.close();
}
