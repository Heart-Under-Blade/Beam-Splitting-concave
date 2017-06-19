#pragma once

#include "geometry_lib.h"
#include "PhysMtr.hpp"
#include "Tracing.h"
#include "CalcTimer.h"
#include "Mueller.hpp"

#define MAX_GROUP_NUM 1024

struct Contribution
{
	Contribution()
		: scatMatrix(0, 0, 0, 0), back(4, 4), forw(4, 4)
	{
		scatMatrix = Arr2D(1, 180 + 1/*TODO: this is 'thetaNum',
				   should i do smth with it?*/, 4, 4);
		scatMatrix.ClearArr();
	}

	Arr2D scatMatrix;	///< Scattering matrix
	matrix back;		///< Mueller matrix in backward direction
	matrix forw;		///< Mueller matrix in forward direction
};

struct TrackGroup
{
	int groupID;
	long long int arr[MAX_GROUP_NUM];
	int size = 0;
	std::vector<std::vector<int>> tracks;

	std::string CreateGroupName() const
	{
		std::string subname;

		if (size <= 2)
		{
			for (int i = 0; i < size; ++i)
			{
				for (int index : tracks[i])
				{
					subname += std::to_string(index) + '_';
				}

				subname += '_' + std::to_string(groupID);
			}
		}

		subname += "_gr_" + std::to_string(groupID);
		return subname;
	}
};

class Tracks : public std::vector<TrackGroup>
{
public:
	int GetMaxGroupID() const // REF: вызывать в конструкторе один раз
	{
		int maxGroupID = 0;

		for (size_t i = 0; i < size(); ++i)
		{
			if ((*this)[i].groupID > maxGroupID)
			{
				maxGroupID = (*this)[i].groupID;
			}
		}

		return ++maxGroupID;
	}

	int FindGroup(long long int trackID) const
	{
		for (size_t i = 0; i < size(); ++i)
		{
			for (int j = 0; j < (*this)[i].size; ++j)
			{
				if ((*this)[i].arr[j] == trackID)
				{
					return (*this)[i].groupID;
				}
			}
		}

		if (size() == 0)
		{
			return 0;
		}

		return -1;
	}

	static void RecoverTrack(const Beam &beam, int facetNum,
							 std::vector<int> &track)
	{
		int coef = facetNum + 1;
		std::vector<int> tmp_track;

		int tmpId = beam.id/coef;
		for (int i = 0; i <= beam.level; ++i)
		{
			int tmp = tmpId%coef;
			tmpId -= tmp;
			tmpId /= coef;
			tmp -= 1;
			tmp_track.push_back(tmp);
		}

		for (int i = tmp_track.size()-1; i >= 0; --i)
		{
			track.push_back(tmp_track.at(i));
		}
	}
};

struct AngleRange
{
	int count;
	double norm;
	double step;

	AngleRange(double begin, double end, int count, double normCoef)
		: count(count)
	{
		norm = normCoef/count;
		step = DegToRad(end - begin)/count;
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

	void TraceIntervalGO(int betaNumber, int gammaNumber, const Tracks &tracks);
	void TraceIntervalGO(int betaNumber, int gammaNumber);
	void TraceSingleOrGO(const double &beta, const double &gamma,
						 const Tracks &tracks);

	void TraceRandomPO(int betaNumber, int gammaNumber, const Cone &bsCone,
					   const Tracks &tracks, double wave);
	void TraceBackScatterPointPO(int betaNumber, int gammaNumber,
								 const Tracks &tracks, double wave);
	void TraceIntervalPO2(int betaNumber, int gammaNumber, const Cone &bsCone,
						  const Tracks &tracks, double wave);
	void TraceSingleOrPO(const double &beta, const double &gamma,
						 const Cone &bsCone, const Tracks &tracks, double wave);

	void setIsCalcOther(bool value); // REF: заменить

private:
	Tracing *m_tracing;

	Symmetry m_symmetry;

	// result scattering martices
	Contribution m_totalMtrx;
	std::vector<Contribution> m_sepatateMatrices; // матрицы для вклада в отдельные траектории

	Point3f m_incidentDir;
	Point3f m_polarizationBasis;
	std::string m_resultDirName;
	std::vector<Arr2DC> J; // Jones matrices

	// light energy balance
	double m_incomingEnergy;
	double m_outcomingEnergy;
	time_t m_startTime;

	// REF: заменить
	bool isCalcOther = false;
	double gNorm;
	Arr2D Other;
	Arr2D All;

	std::string m_statistics;

private:
	void CleanJ(int maxGroupID, const Cone &bsCone);
	void HandleBeamsGO(std::vector<Beam> &outBeams, double beta);
	void HandleBeamsGO(std::vector<Beam> &outBeams, double beta, const Tracks &tracks);
	void HandleBeamsPO(std::vector<Beam> &outBeams, const Cone &bsCone, double wavelength, const Tracks &tracks);
	void HandleBeamsPO2(std::vector<Beam> &outBeams, const Cone &bsCone, double wavelength, int groupID);
	void HandleBeamsBackScatterPO(std::vector<Beam> &outBeams,
								  double wavelength, const Tracks &tracks);
	void SetJnRot(Beam &beam, const Point3f &T,
				  const Point3d &vf, const Point3d &vr, matrixC &Jn_rot);
	void AddResultToMatrix(Arr2D &M, int maxGroupID, const Cone &bsCone,
							  double norm = 1);
	void AddResultToMatrices(std::vector<Arr2D> &M, int maxGroupID,
									 const Cone &bsCone, double norm = 1);
	void WriteConusMatrices(std::ofstream &outFile, const Arr2D &sum,
						const Cone &bsCone);
	void AddToSumMatrix(const Cone &bsCone, double norm, int q, Arr2D &M_);
	void OutputProgress(int betaNumber, long long count, CalcTimer &timer);
	void ExtractPeaksGO(int EDF, double NRM, Contribution &contr);
	void WriteResultsToFileGO(double NRM, const std::string &filename,
							  Contribution &contr);
	std::string GetFileName(const std::string &filename);

	double CalcNorm(long long orNum);
	double CalcTotalScatteringEnergy();
	void RotateMuller(const Point3f &dir, matrix &bf);
	void AddToResultMullerGO(const Point3f &dir, matrix &bf, double area,
							 Contribution &contr);
	void WriteResultToSeparateFilesGO(double NRM, int EDF, const Tracks &tracks);
	void AllocGroupMatrices(std::vector<Arr2D> &mtrcs, size_t maxGroupID);

	void CreateGroupResultFiles(const Tracks &tracks, const std::string &dirName,
								std::vector<std::ofstream*> &groupFiles);
	void AllocJ(int m, int n, int size);
	void CleanJ();
	void OutputStartTime(CalcTimer timer);
	void OutputStatisticsGO(int orNumber, double D_tot, double NRM,
						  CalcTimer &timer);
};
