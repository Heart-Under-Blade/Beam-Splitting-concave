#include "TracerPOTotal.h"

#include "HandlerPO.h"
#include "Mueller.hpp"
#include "PhysMtr.hpp"
#include "MullerMatrix.h"
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
	long long nOrientations = (betaRange.number + 1) * gammaRange.number;
	ofstream outFile(m_resultDirName + "_out.dat", ios::out);

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
//	for (int i = 6; i <= 6/*betaRange.number*/; ++i)
	{
		beta = i*betaRange.step;

//#ifdef _DEBUG // DEB
//		beta = DegToRad(170);
//#endif
		double sinZenith = (i == 0 || i == betaRange.number)
				? (1.0-cos(0.5*betaRange.step))/normIndex
				: (cos((i-0.5)*betaRange.step) -
				   cos((i+0.5)*betaRange.step))/normIndex;

//		std::cout << sinZenith << std::endl;
		m_handler->SetSinZenith(sinZenith);

		for (int j = 0; j < gammaRange.number; ++j)
//		for (int j = 65; j < 67/*gammaRange.number*/; ++j)
		{
			gamma = j*gammaRange.step;

			m_particle->Rotate(M_PI-beta, M_PI+gamma, 0);
			m_scattering->FormShadowBeam(outBeams);
			m_scattering->ScatterLight(M_PI-beta, M_PI+gamma, outBeams);

			m_handler->HandleBeams(outBeams);

			outBeams.clear();

			OutputProgress(nOrientations, count, i, j, timer);
			++count;

#ifdef _DEBUG // DEB
//			if (i == 9)
//			{
				double m = static_cast<HandlerPO*>(m_handler)->M(0,100)[0][0];
				m_handler->m_logFile << i << ' ' << j << ' ' << m << endl;
//			}
#endif
		}
	}

	std::string dir = CreateFolder(m_resultDirName);
	m_resultDirName = dir + m_resultDirName + '\\' + m_resultDirName;
	m_handler->WriteMatricesToFile(m_resultDirName);

//	std::ofstream outFile1(m_resultDirName + ".dat", std::ios::out);
//	for (int i = 0; i <= m_handler->m_sphere.nAzimuth; ++i)
//	{
//		for (int j = 0; j <= m_handler->m_sphere.nZenith; ++j)
//		{
//			matrix m(4,4);
//			m.Fill(0);
//			m = static_cast<HandlerPO*>(m_handler)->M(i,j);
//			outFile1 << i << " " << j << " ";
//			outFile1 << m << endl;
//		}
//	}

	OutputLogPO(timer, nOrientations, dir);
	outFile.close();
}
