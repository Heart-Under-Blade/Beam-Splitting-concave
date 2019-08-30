#include "TracerPO.h"

using namespace std;

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

	CalcTimer timer;
	OutputStartTime(timer);

	ofstream outFile(m_resultDirName, ios::out);

	int halfGammaNum = range.nAzimuth/2;

	for (int i = 0; i <= range.nZenith; ++i)
	{
		angle.zenith = i*range.step.zenith;

		for (int j = -halfGammaNum; j <= halfGammaNum; ++j)
		{
			angle.azimuth = j*range.step.azimuth;

			m_particle->Rotate(angle);
			m_scattering->ScatterLight(outBeams);
			m_handler->HandleBeams(outBeams);
			outBeams.clear();
		}

		m_handler->WriteMatricesToFile(m_resultDirName);

		OutputProgress(nOrientations, i, timer);
	}

	outFile.close();
}
