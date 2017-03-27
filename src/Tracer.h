#pragma once

#include "geometry_lib.h"
#include "PhysMtr.hpp"
#include "Tracing.h"

struct TrackGroup
{
	int groupID;
	long long int arr[1024];
	int size = 0;
};

struct Tracks
{
	TrackGroup groups[32];
	int count = 0;

	int GetMaxGroupID()
	{
		int maxGroupID = 0;

		for (int i = 0; i < count; ++i)
		{
			if (groups[i].groupID > maxGroupID)
			{
				maxGroupID = groups[i].groupID;
			}
		}

		return maxGroupID;
	}

	int GetGroupID(long long int trackID)
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

	double GetStep()
	{
		return DegToRad(end - begin)/count;
	}
};

struct Cone
{
	Cone(double radius, int phiCount, int thetaCount)
		: radius(radius), phiCount(phiCount), thetaCount(thetaCount)
	{
		dPhi = M_2PI/phiCount;
		dTheta = radius/thetaCount;
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
	Tracer(const Tracing *tracing, const std::string resultFileName);

	void TraceIntervalPO(const AngleInterval &betaI, const AngleInterval &gammaI,
						 const Cone &bsCone, const Tracks &tracks, double wave);
	void TraceIntervalGO(const AngleInterval &betaI, const AngleInterval &gammaI);
	void TraceSingleOr(const double &beta, const double &gamma);

private:
	std::string m_resultFileName;
	Point3f m_polarizationBasis;
	Point3f m_incidentDir;
	Tracing *m_tracing;
	Arr2D m_mxd;
	double m_gammaNorm;
	std::vector<Arr2DC> J; // Jones matrices

private:
	void CleanJ(int maxGroupID, const Cone &bsCone);
	void HandleBeamsPO(std::vector<Beam> &outBeams, const Cone &bsCone, double wavelength, const Tracks &tracks);
	void SetJnRot(Beam &beam, const Point3f &T,
				  const Point3d &vf, const Point3d &vr, matrixC &Jn_rot);
	void AddResultToSumMatrix(Arr2D &M_, int maxGroupID, const Cone &bsCone,
							  const AngleInterval &gammaI);
	void WriteSumMatrix(std::ofstream &outFile, const Arr2D &sum);
};
