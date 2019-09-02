#include "TracerPO.h"

using namespace std;

<<<<<<< HEAD
TracerPO::TracerPO(Particle *particle, int nActs, const string &resultFileName)
	: Tracer(particle, nActs, resultFileName)
{
}

void TracerPO::TraceRandom(const AngleRange &betaRange, const AngleRange &gammaRange)
{
=======
TracerPO::TracerPO(Particle *particle, Scattering *scattering,
				   const string &resultFileName)
	: Tracer(particle, scattering, resultFileName)
{
}

void TracerPO::TraceRandom(const OrientationRange &range)
{
	vector<Beam> outBeams;
	Orientation angle;
	long long nOrientations = range.nZenith * range.nAzimuth;

>>>>>>> origin/refactor
	CalcTimer timer;
	long long count = 0;

	ofstream outFile(m_resultDirName, ios::out);

<<<<<<< HEAD
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
=======
	int halfGammaNum = range.nAzimuth/2;

	for (int i = 0; i <= range.nZenith; ++i)
	{
		angle.zenith = i*range.step.zenith;

		for (int j = -halfGammaNum; j <= halfGammaNum; ++j)
		{
			angle.azimuth = j*range.step.azimuth;
>>>>>>> origin/refactor

			m_particle->Rotate(beta, gamma, 0);
			m_scattering->ScatterLight(outBeams);

			m_handler->HandleBeams(outBeams);
			outBeams.clear();
		}

		m_handler->WriteMatricesToFile(m_resultDirName);

<<<<<<< HEAD
		OutputProgress(betaRange.number, count, timer);
		++count;
=======
		OutputProgress(nOrientations, i, timer);
>>>>>>> origin/refactor
	}

	outFile.close();
}

void TracerPO::TraceFixed(const double &beta, const double &gamma)
{
	ofstream outFile(m_resultDirName, ios::out);
	vector<Beam> outBeams;

	double b = DegToRad(beta);
	double g = DegToRad(gamma);
	m_particle->Rotate(b, g, 0);
	m_scattering->ScatterLight(outBeams);

	m_handler->HandleBeams(outBeams);
	outBeams.clear();
	m_handler->WriteMatricesToFile(m_resultDirName);
	outFile.close();
}
