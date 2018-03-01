#pragma once

#include "Beam.h"
#include "Scattering.h"
#include "PhysMtr.hpp"

struct ContributionGO
{
	ContributionGO()
		: scatMatrix(0, 0, 0, 0), back(4, 4), forw(4, 4)
	{
		scatMatrix = Arr2D(1, 180 + 1/*TODO: this is 'thetaNum',
				   should i do smth with it?*/, 4, 4);
		scatMatrix.ClearArr();
		back.Fill(0);
		forw.Fill(0);
	}

	Arr2D scatMatrix;	///< Scattering matrix
	matrix back;		///< Mueller matrix in backward direction
	matrix forw;		///< Mueller matrix in forward direction
};

class HandlerGO
{
public:
	HandlerGO(Particle *particle, Light *incidentLight, float wavelength = 0);

	virtual void HandleBeams(std::vector<Beam> &beams, double angle);
	virtual void WriteMatricesToFile(double norm, std::string &filename);
	void SetTracks(Tracks *tracks);
	double CalcTotalScatteringEnergy(double norm);
	void SetAbsorbtionAccounting(bool value);

protected:
	Particle *m_particle;
	Light *m_incidentLight;
	Tracks *m_tracks;

	bool m_isAccountAbsorbtion;
	float m_wavelength;

	ContributionGO m_totalContrib;	// result scattering martices contribution
	std::vector<ContributionGO> m_groupContrib; // group contibution of beams

	std::ofstream m_logFile;

protected:
	double BeamCrossSection(const Beam &beam) const;
	void AddToResultMullerGO(const Point3f &dir, matrix &bf, double area,
							 ContributionGO &contr);
	void RotateMuller(const Point3f &dir, matrix &bf);
	void AverageOverAlpha(int EDF, double norm, ContributionGO &contrib);

	void WriteToFile(ContributionGO &contrib, double norm,
					 const std::string &filename);
	double CalcOpticalPathAbsorption(const Beam &beam);
	Point3f CalcK(const Beam &beam, std::vector<int> &tr);

private:
	void ExtractPeaks(double *b, double *f, double norm);
};


class HandlerTotalGO : public HandlerGO
{
public:
	HandlerTotalGO(Particle *particle, Light *incidentLight, float wavelength = 0);

	void HandleBeams(std::vector<Beam> &beams, double angle) override;
	void WriteMatricesToFile(double norm, std::string &filename) override;
};


class HandlerGroupGO : public HandlerGO
{
public:
	HandlerGroupGO(Particle *particle, Light *incidentLight, float wavelength = 0);

	void HandleBeams(std::vector<Beam> &beams, double angle) override;
	void WriteMatricesToFile(double norm, std::string &dir) override;
};
