#pragma once

#include "geometry_lib.h"
#include "PhysMtr.hpp"
#include "Tracing.h"

struct TrackGroup
{
	int groupID;
	long long int arr[1024];
	int size = 0;
	std::vector<std::vector<int>> tracks;
};

struct Tracks
{
	TrackGroup groups[32];
	int count = 0;

	int GetMaxGroupID() const // REF: вызывать в конструкторе
	{
		int maxGroupID = 0;

		for (int i = 0; i < count; ++i)
		{
			if (groups[i].groupID > maxGroupID)
			{
				maxGroupID = groups[i].groupID;
			}
		}

		return ++maxGroupID;
	}

	int GetGroupID(long long int trackID) const
	{
		for (int i = 0; i < count; ++i)
		{
			for (int j = 0; j < groups[i].size; ++j)
			{
				if (groups[i].arr[j] == trackID)
				{
					return groups[i].groupID;
				}
			}
		}

		return -1;
	}
};

struct AngleInterval
{
	double begin;
	double end;
	int count;
	double norm;

	void SetNorm(const double &coef)
	{
		norm = coef/count;
	}

	double GetStep() const
	{
		return DegToRad(end - begin)/count;
	}
};

/**
 * @brief The Cone struct
 * Backscattering cone divided by cells
 */
struct Cone
{
	Cone(double radius, int phiCount, int thetaCount)
		: radius(radius), phiCount(phiCount), thetaCount(thetaCount)
	{
		dPhi = M_2PI/(phiCount+1);
		dTheta = DegToRad(radius)/thetaCount;
	}

	double radius;
	int phiCount;
	int thetaCount;
	double dPhi;
	double dTheta;
};

class Tracer
{
public:
	Tracer(Tracing *tracing, const std::string resultFileName);

	void TraceIntervalPO(const AngleInterval &betaI, const AngleInterval &gammaI,
						 const Cone &bsCone, const Tracks &tracks, double wave);
	void TraceIntervalPO2(const AngleInterval &betaI, const AngleInterval &gammaI,
						 const Cone &bsCone, const Tracks &tracks, double wave);
	void TraceIntervalGO(const AngleInterval &betaI, const AngleInterval &gammaI);
	void TraceSingleOrPO(const double &beta, const double &gamma,
						 const Cone &bsCone, const Tracks &tracks, double wave);

private:
	Tracing *m_tracing;
	Arr2D m_mxd;
	Point3f m_incidentDir;
	Point3f m_polarizationBasis;
	std::string m_resultFileName;
	double m_gammaNorm;
	std::vector<Arr2DC> J; // Jones matrices

private:
	void CleanJ(int maxGroupID, const Cone &bsCone);
	void HandleBeamsPO(std::vector<Beam> &outBeams, const Cone &bsCone, double wavelength, const Tracks &tracks);
	void HandleBeamsPO2(std::vector<Beam> &outBeams, const Cone &bsCone, double wavelength, int groupID);
	void SetJnRot(Beam &beam, const Point3f &T,
				  const Point3d &vf, const Point3d &vr, matrixC &Jn_rot);
	void AddResultToSumMatrix(Arr2D &M_, int maxGroupID, const Cone &bsCone,
							  double norm = 1);
	void WriteSumMatrix(std::ofstream &outFile, const Arr2D &sum,
						const Cone &bsCone);
	void EraseConsoleLine(int lenght); // REF: вынести куда-нибудь
	void AddToSumMatrix(const Cone &bsCone, double norm, int q, Arr2D &M_);
	void PrintProgress(int betaNumber, long long count);
};
