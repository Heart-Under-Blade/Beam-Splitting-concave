#include "TracerPO.h"

using namespace std;

TracerPO::TracerPO(Particle *particle, Scattering *scattering,
				   const string &resultFileName)
	: LightTracer(particle, scattering, resultFileName)
{
}

void TracerPO::TraceRandom(const AngleRange &zenithRange,
						   const AngleRange &azimuthRange)
{
	vector<Beam> outBeams;
	Orientation angle;

	CalcTimer timer;
	OutputStartTime(timer);

	ofstream outFile(m_resultDirName, ios::out);

	int halfGammaNum = azimuthRange.number/2;

	for (int i = 0; i <= zenithRange.number; ++i)
	{
		angle.zenith = i*zenithRange.step;

		for (int j = -halfGammaNum; j <= halfGammaNum; ++j)
		{
			angle.azimuth = j*azimuthRange.step;

			m_particle->Rotate(angle);
			m_scattering->ScatterLight(outBeams);
			m_handler->HandleBeams(outBeams);
			outBeams.clear();
		}

		m_handler->WriteMatricesToFile(m_resultDirName);

		OutputProgress(zenithRange.number, i, timer);
	}

	outFile.close();
}
