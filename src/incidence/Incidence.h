#pragma once

#include "PathedBeam.h"
#include "Splitting.h"

template <class T>
class SplittedBeams
{
public:
	T internal;
	T external;
};

class Incidence
{
public:
	Incidence(const complex &ri);
	virtual void ComputeDirections(Beam &/*beam*/, SplittedBeams<Beam> &/*beams*/) const {};
	virtual void ComputeJonesMatrices(Beam &/*beam*/, SplittedBeams<Beam> &/*beams*/) const {};

	virtual void ComputeOpticalPaths(const PathedBeam &beam,
									 SplittedBeams<PathedBeam> &beams) const;

	void ComputeReflectedDirection(Vector3f &dir) const
	{
		dir = r - m_normal;
		Point3f::Normalize(dir); // REF, OPT: нужно ли это нормализовать всегда?
	}

	void ComputeRefractedDirection(Vector3f &dir) const
	{
		dir = r/sqrt(s) + m_normal;
		Point3f::Normalize(dir); // REF, OPT: нужно ли это нормализовать всегда?
	}

public:
	double reRiEff;
	double cosA;
	Point3f m_normal;
	double s;
	Point3f r;

protected:
	complex m_ri;
};
