#pragma once

#include "Beam.h"
#include "Particle.h"
#include "Intersection.h"

#include <vector>

#define MAX_BEAM_DEPT 128

struct BeamInfo
{
	Beam beam;
	int lastFacetIndex;
	int dept;
};

struct BeamStack
{
public:
	BeamStack()
	{
		m_currentDept = 0;
	}

	void Push(const Beam &beam, int lastFacetIndex, int beamDept)
	{
		m_beams[m_currentDept] = beam;
		m_lastFacetIndices[m_currentDept] = lastFacetIndex;
		m_beamDepts[m_currentDept] = beamDept;

		++m_currentDept;
	}

	void Pop(Beam &beam, int &lastFacetIndex, int &beamDept)
	{
		--m_currentDept;

		beam = m_beams[m_currentDept];
		lastFacetIndex = m_lastFacetIndices[m_currentDept];
		beamDept = m_beamDepts[m_currentDept];
	}

	bool IsEmpty()
	{
		return m_currentDept == 0;
	}

private:
	int m_currentDept;

	Beam m_beams[MAX_BEAM_DEPT];
	int m_lastFacetIndices[MAX_BEAM_DEPT];
	int m_beamDepts[MAX_BEAM_DEPT];
};

struct OutBeam
{
	Beam beam;

	int track[MAX_BEAM_DEPT];
	int trackSize;

	OutBeam(const Beam &p_beam, const int *p_track, int p_trackSize)
	{
		beam = p_beam;
		trackSize = p_trackSize;

		for (int i = 0; i < trackSize; ++i)
		{
			track[i] = p_track[i];
		}
	}
};

class Tracing
{
public:
	Tracing(Particle *p_particle,
			bool p_isOpticalPath,
			const Point3f &p_incidentDir,
			const Point3f &p_polarizationBasis,
			int p_interReflectionNumber);

	double TraceParticle(double beta, double gamma);
	void ClearTracing();

public:
	std::vector<OutBeam> outcomingBeams;

private:

	/**
	 * @brief The IncidenceCase enum
	 * Case of beam incidence to facet of particle
	 */
	enum class IncidenceCase : int
	{
		Normal,
		Slopping
	};

	Particle *particle;			///> scattering particle (crystal)
	Point3f incidentDirection;	///> first incedent beam direction
	Point3f polarizationBasis;	///>
	bool isOpticalPath;
	int interReflectionNumber;

	int track[MAX_BEAM_DEPT]; ///> path of the current beam (numbers of facets refracted the beam)
	int trackSize = 0;

	const double FAR_ZONE_DISTANCE = 10000.0;
	const double LOW_ENERGY_LEVEL = 2e-12;

private:
	void TraceInternal(BeamInfo *tree, int treeDept);

	void SetBeam(Beam &beam, const Beam &other, const Point3f &dir, const Point3f &e,
				 const complex &coef1, const complex &coef2) const;

	bool Intersect(int facetIndex, const Beam& originBeam, Beam &intersectBeam) const;

	void SplitIncidentDirection(const Point3f &incidentDirection, double cosIncident, int normalIndex,
								Point3f &reflectionDir, Point3f &refractionDir) const;

	void SetBeamShapesByFacet(int facetIndex, Beam &inBeam, Beam &outBeam) const;

	void SetOutputBeam(__m128 *_output_points, int outputSize, Beam &outputBeam) const;

	bool ProjectToFacetPlane(const Beam& inputBeam, __m128 *_output_points,
					__m128 _normal, int facetIndex) const;

	void CalcOpticalPathInternal(double Nr, const Beam &incidentBeam, Beam &outBeam, Beam &inBeam) const;

	void SplitBeam(const Beam &incidentBeam, Beam &inBeam, Beam &outBeam, double Nr,
				   IncidenceCase incidenceCase);

	void InvertBeamPointOrder(Beam &outBeam, const Beam &inBeam);

	inline bool isEnough(const BeamInfo &info);
	inline void changeTrack(int &lastBeamDept, const BeamInfo &info);
};
