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

	double normIndex = range.step.azimuth/(range.to.azimuth - range.from.azimuth);
	m_handler->SetNormIndex(normIndex);

	timer.Start();

	for (int i = 0; i <= range.nZenith; ++i)
	{
		angle.zenith = i*range.step.zenith;

		for (int j = -halfGammaNum; j <= halfGammaNum; ++j)
		{
			angle.azimuth = j*range.step.azimuth;

			m_particle->Rotate(angle);
			m_scattering->Scatter(&outBeams);
			m_handler->HandleBeams(outBeams);
			outBeams.clear();
		}

		m_handler->WriteMatricesToFile(m_resultDirName);

		OutputProgress(nOrientations, i, timer);
	}

	outFile.close();
}

void TracerPO::TraceFixed(const Orientation &orientation)
{
	ofstream outFile(m_resultDirName, ios::out);
	vector<Beam> outBeams;

	m_particle->Rotate(orientation.ToRadians());
	m_scattering->ExtractShadowBeam();
	m_scattering->Scatter(&outBeams);

	m_handler->HandleBeams(outBeams);
	outBeams.clear();

	m_handler->WriteMatricesToFile(m_resultDirName);
	outFile.close();
}
