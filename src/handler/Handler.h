#pragma once

#include "Beam.h"
#include "Scattering.h"
#include "PhysMtr.hpp"
#include "MullerMatrix.h"
#include "Tracks.h"
#include "ScatteringFiles.h"

#define SPHERE_RING_NUM		180		// number of rings no the scattering sphere

/**
 * @brief The Cone struct
 * Backscattering cone divided by cells
 */
struct Conus
{
	Conus(double radius, int phiCount, int thetaCount)
		: radius(radius), phiCount(phiCount), thetaCount(thetaCount)
	{
		dPhi = M_2PI/(phiCount+1);
		dTheta = Orientation::DegToRad(radius)/thetaCount;
	}

	double radius;
	int phiCount;
	int thetaCount;
	double dPhi;
	double dTheta;
};


class PointContribution
{
public:
	PointContribution(int nGroups, double normIndex)
		: m_nGroups(nGroups),
		  m_normIndex(normIndex)
	{
		groupJones.resize(nGroups);
		groupMuellers.resize(nGroups);
		ResetJones();
	}

	void AddToMueller(const Matrix2x2c &jones)
	{
		MuellerMatrix m(jones);
		m *= m_normIndex;
		rest += m;
	}

	void AddToGroup(const Matrix2x2c &jones, int groupId)
	{
		groupJones[groupId] += jones;
	}

	void SumGroupTotal()
	{
		for (int gr = 0; gr < m_nGroups; ++gr)
		{
			MuellerMatrix m(groupJones[gr]);
			m *= m_normIndex;
			groupMuellers[gr] += m;
			groupTotal += m;
		}

		ResetJones();
	}

	void SumTotal()
	{
		total += groupTotal;
		total += rest;
	}

	const MuellerMatrix &GetGroupTotal() const
	{
		return groupTotal;
	}

	const MuellerMatrix &GetTotal() const
	{
		return total;
	}

	const MuellerMatrix &GetRest() const
	{
		return rest;
	}

	const MuellerMatrix &GetGroupMueller(int groupID)
	{
		return groupMuellers.at(groupID);
	}

	void Reset()
	{
		groupTotal.Reset();
		rest.Reset();
		total.Reset();

		for (MuellerMatrix &m : groupMuellers)
		{
			m.Reset();
		}
	}

private:
	std::vector<Matrix2x2c> groupJones;
	std::vector<MuellerMatrix> groupMuellers;

	MuellerMatrix groupTotal;
	MuellerMatrix rest;
	MuellerMatrix total;

	int m_nGroups;
	double m_normIndex;

	void ResetJones()
	{
		for (Matrix2x2c &j : groupJones)
		{
			j.Fill(0.f);
		}
	}
};

class Handler // REF: rename to "Detector"
{
public:
	Handler(Particle *particle, Light *incidentLight, float wavelength = 0);

	virtual void HandleBeams(std::vector<Beam> &beams);
	virtual void SetTracks(Tracks *tracks);
	void SetScattering(Scattering *scattering);
	virtual void WriteMatricesToFile(std::string &destName);

	void SetAbsorbtionAccounting(bool value);
	void SetNormIndex(double normIndex);

	Light *m_incidentLight;
	Tracks *m_tracks;
	float m_wavelength;

	static double BeamCrossSection(const Beam &beam);

protected:
	void ApplyAbsorbtion(Beam &beam);
	void OutputPaths(const Beam &beam, const OpticalPath &path);
	OpticalPath ComputeOpticalPath(const Beam &beam);

protected:
	Scattering *m_scattering;

	Particle *m_particle;
	bool m_hasAbsorbtion;
	double m_normIndex;
	std::ofstream m_logFile;
	std::ofstream m_absLogFile;
	double m_cAbs;
	long long count = 0;
};


