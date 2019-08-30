#include "HandlerPOTotal.h"

#include "Mueller.hpp"
#include <iostream>

HandlerPOTotal::HandlerPOTotal(Particle *particle, Light *incidentLight,
							   double wavelength)
	: HandlerPO(particle, incidentLight, wavelength)
{
	m_Lp = new matrix(4, 4);

	(*m_Lp)[0][0] = 1;
	(*m_Lp)[0][1] = 0;
	(*m_Lp)[0][2] = 0;
	(*m_Lp)[0][3] = 0;

	(*m_Lp)[1][0] = 0;
	(*m_Lp)[1][3] = 0;

	(*m_Lp)[2][0] = 0;
	(*m_Lp)[2][3] = 0;

	(*m_Lp)[3][0] = 0;
	(*m_Lp)[3][1] = 0;
	(*m_Lp)[3][2] = 0;
	(*m_Lp)[3][3] = 1;

	m_Ln = new matrix(4, 4);

	(*m_Ln) = (*m_Lp);
}

void HandlerPOTotal::WriteMatricesToFile(std::string &destName)
{
	std::ofstream outFile(destName + ".dat", std::ios::out);
	matrix sum(4, 4);

	auto &Lp = *m_Lp;
	auto &Ln = *m_Ln;

	int &nT = m_sphere.nZenith;
	double &dT = m_sphere.zenithStep;

	int &nP = m_sphere.nAzimuth;

	for (int t = nT-1; t >= 0; --t)
	{
		sum.Fill(0.0);
		double tt = 180.0 - RadToDeg(t*dT);

		for (int p = 0; p <= nP; ++p)
		{
			double radPhi = -p*m_sphere.azinuthStep;
			matrix m = M(p, t);
//			if (t == 100) {
//				std::cout << p << ' ' << m[0][0] << ' ' << std::endl;
//			}

			Lp[1][1] = cos(2*radPhi);
			Lp[1][2] = sin(2*radPhi);
			Lp[2][1] = -Lp[1][2];
			Lp[2][2] = Lp[1][1];

			Ln[1][2] = -Lp[1][2];
			Ln[2][1] = -Lp[2][1];

			if (t == 0)
			{
				sum += Lp*m*Lp;
			}
			else if (t == nT-1)
			{
				sum += Ln*m*Lp; // OPT: вынести Ln в отдельный случай
			}
			else
			{
				sum += m*Lp;
			}
		}

		double dS2 = (t == 0 || t == (nT-1)) ? 1.0-cos(0.5*dT)
											 : cos((t-0.5)*dT)-cos((t+0.5)*dT);
		dS2 *= M_2PI;
		outFile << std::endl << tt << ' ' << dS2 << ' ';
		outFile << sum/m_sphere.nAzimuth;
	}
}

void HandlerPOTotal::AddToMueller()
{
	for (size_t q = 0; q < m_diffractedMatrices.size(); ++q)
	{
		auto &diffM = m_diffractedMatrices[q];

		for (int t = 0; t < m_sphere.nZenith; ++t)
		{
			for (int p = 0; p < m_sphere.nAzimuth; ++p)
			{
				matrix m = Mueller(diffM(p, t));
#ifndef _DEBUG // DEB
				complex ddd[4];
				ddd[0] = diffM(p, t)[0][0];
				ddd[1] = diffM(p, t)[0][1];
				ddd[2] = diffM(p, t)[1][0];
				ddd[3] = diffM(p, t)[1][1];

				double fff[4];
				fff[0] = m[0][0];
				fff[1] = m[0][1];
				fff[2] = m[1][0];
				fff[3] = m[1][1];
#endif
				m *= m_sinZenith;

//				if (t == 100)
//					m_logFile << m[0][0] << std::endl;

				M.insert(p, t, m);
			}
		}
	}
}
