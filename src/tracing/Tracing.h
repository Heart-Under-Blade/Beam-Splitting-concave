#pragma once

#include "Beam.h"
#include "Particle.h"
#include "Intersection.h"

#include <float.h>
#include <vector>

//#define MAX_BEAM_REFL_NUM 32768
#define MAX_BEAM_REFL_NUM 65536
//#define MAX_BEAM_REFL_NUM 1048576

#define EPS_M_COS_90	-1.7453292519943295769148298069306e-10	//cos(89.99999999)
#define EPS_COS_90		1.7453292519943295769148298069306e-10	//cos(89.99999999)
#define EPS_COS_00		0.99999999998254670756866631966593		//1- cos(89.99999999)
#define MAX_GROUP_NUM	1024

/**
 * @brief The BeamTree struct
 * Tree of beams (works like stack).
 */
struct BeamTree
{
	Beam tree[MAX_BEAM_REFL_NUM];
	int size = 0;

//	void Push(const Beam &beam)
//	{
//		tree[size++] = beam;
//	}
};

struct TrackGroup
{
	int groupID;
	BigInteger arr[MAX_GROUP_NUM];
	int size = 0;
	std::vector<std::vector<int>> tracks;

	std::string CreateGroupName() const
	{
		std::string subname;
		subname += "gr_" + std::to_string(groupID);
		return subname;
	}
};

class Tracks : public std::vector<TrackGroup>
{
public:
	int FindGroup(const BigInteger &trackID) const
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

		BigInteger tmpId = beam.id/coef;
		for (int i = 0; i <= beam.level; ++i)
		{
			int tmp = (tmpId%coef).toInt();
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

class Tracing // REF: поменять название на ~Splitter
{
public:
	// REF: убрать в протектед потом
	Particle *m_particle;			///< scattering particle (crystal)
protected:
	Facet *m_facets;

	Light *m_incidentLight;
	Point3f m_incidentDir;
	Point3f m_polarBasis;

	bool m_isOpticalPath;
	int m_interReflectionNumber;

//	std::vector<Beam> m_beamTree;
	Beam m_beamTree[MAX_BEAM_REFL_NUM];	///< tree of beams (works like stack)
	int m_treeSize;
	double m_incommingEnergy;

	const double FAR_ZONE_DISTANCE = 10000.0;
	const double LOW_ENERGY_LEVEL = 2e-12;

private:
	double ri_coef_re;
	double ri_coef_im;
	complex m_refrIndex;

public:
	Tracing(Particle *particle, Light *incidentLight, bool isOpticalPath,
			int interReflectionNumber);

	virtual void SplitBeamByParticle(double /*beta*/, double /*gamma*/, std::vector<Beam> &/*scaterredBeams*/) {}
	virtual void SplitBeamByParticle(double beta, double gamma, const std::vector<std::vector<int>> &tracks,
									 std::vector<Beam> &scaterredBeams);

	double GetIncomingEnergy() const;

//	double CrossSection(const Point3f &beamDir) const;

protected:
	void SetBeamOpticalParams(int facetId, Beam &inBeam, Beam &outBeam);

	void SetBeam(Beam &beam, const Beam &other, const Point3f &dir,
				 const complex &coef1, const complex &coef2) const;

	void Difference(const Polygon &subject, const Point3f &subjNormal,
					const Polygon &clip, const Point3f &clipNormal,
					const Point3f &clipDir, Polygon *difference, int &resultSize) const;

	bool Intersect(int facetId, const Beam& beam, Polygon &intersection) const;

	void SetPolygonByFacet(int facetId, Polygon &polygon) const;

	void CalcOpticalPath(double cosIN, const Beam &incidentBeam,
								 Beam &inBeam, Beam &outBeam) const;

	bool IsTerminalBeam(const Beam &beam);

	void TraceFirstBeam(int facetId, Beam &inBeam, Beam &outBeam);

	void TraceSecondaryBeams(Beam &incidentBeam, int facetID,
								  Beam &inBeam, std::vector<Beam> &outBeams);

	void RotatePolarisationPlane(const Point3f &dir, const Point3f &facetNormal,
								 Beam &beam);

	void SetSloppingBeamParams_initial(const Point3f &beamDir, double cosIN, int facetId,
									   Beam &inBeam, Beam &outBeam);

	void SetNormalIncidenceBeamParams(double cosIN, const Beam &incidentBeam,
									  Beam &inBeam, Beam &outBeam);

	void SetSloppingIncidenceBeamParams(double cosIN, const Point3f &normal,
										Beam &incidentBeam, Beam &inBeam, Beam &outBeam,
										bool &isTrivialIncidence);

	void CalcFacetEnergy(int facetID, const Polygon &lightedPolygon);

	void CalcOpticalPath_initial(Beam &inBeam, Beam &outBeam);

	void PushBeamToTree(Beam &beam, int facetId, int level, Location location);
	void PushBeamToTree(Beam &beam, int facetId, int level);
	void PushBeamToTree(Beam &beam);

	void SetBeamID(Beam &beam);

private:
	double CalcReRI(const double &cosIN) const;

	void SetTrivialIncidenceBeamParams(double cosIN, double Nr, const Point3f &normal,
									   Point3f r0, double s, const Beam &incidentBeam,
									   Beam &inBeam, Beam &outBeam);

	void SetCompleteReflectionBeamParams(double cosIN, double Nr, const Beam &incidentBeam,
										 Beam &inBeam);

	void DivideBeamDirection(const Point3f &incidentDir, double cosIN, const Point3f &normal,
							 Point3f &reflDir, Point3f &refrDir) const;

	void SetOutputPolygon(__m128 *_output_points, int outputSize,
						  Polygon &polygon) const;

	bool ProjectToFacetPlane(const Polygon &polygon, const Point3f &dir,
							 const Point3f &normal, __m128 *_projection) const;
};
