#pragma once

#include "Handler.h"

#define SPHERE_RING_NUM	180 // number of rings no the scattering sphere
#define BIN_SIZE M_PI/SPHERE_RING_NUM

class ContributionGO
{
public:
	ContributionGO() : muellers(0, 0, 0, 0), back(4, 4), forward(4, 4)
	{
		muellers = Arr2D(1, 180 + 1/*TODO: this is 'thetaNum',
				   should i do smth with it?*/, 4, 4);
		muellers.ClearArr();
		back.Fill(0);
		forward.Fill(0);
	}

	void AddMueller(float fAngle, int angle, const matrix &m)
	{
		if (fAngle >= 180-FLT_EPSILON)
		{
			forward += m;
		}
		else if (fAngle <= FLT_EPSILON)
		{
			back += m;
		}
		else
		{
			muellers.insert(0, angle, m);
		}
	}

	Arr2D muellers;		///< Scattering matrices
	matrix back;		///< Mueller matrix in backward direction
	matrix forward;		///< Mueller matrix in forward direction
};


class HandlerGO : public Handler
{
public:
	HandlerGO(Particle *particle, Light *incidentLight, double wavelength = 0);

	void SetTracks(Tracks *tracks) override;

	double ComputeTotalScatteringEnergy();
	void WriteLog(const std::string &str);

	void MultiplyMueller(const Beam &beam, matrix &m);
	virtual void WriteMatricesToFile(std::string &destName) override;

protected:
	ContributionGO m_totalContrib;	// result scattering martices contribution
	std::vector<ContributionGO> m_tracksContrib; // group contibution of beams

protected:
	matrix ComputeMueller(float zenAng, Beam &beam);
	void RotateMuller(const Point3f &dir, matrix &bf);
	void AverageOverAlpha(int EDF, double norm, ContributionGO &contrib,
						  const std::string &destDir);

	void WriteToFile(ContributionGO &contrib, double norm,
					 const std::string &filename);

	double ComputeOpticalPathAbsorption(const Beam &beam);
	Point3f CalcK(std::vector<int> &tr);

private:
	void ExtractPeaks(double *b, double *f, double norm,
					  const std::string &destDir);
};
