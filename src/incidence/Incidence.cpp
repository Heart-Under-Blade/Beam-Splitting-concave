#include "Incidence.h"

Incidence::Incidence(const complex &ri)
	: m_ri(ri)
{
}

void Incidence::ComputeOpticalPaths(const PathedBeam &beam,
									SplittedBeams<PathedBeam> &beams) const
{
	if (beam.opticalPath < FLT_EPSILON)
	{
		Point3f p = beams.internal.Center();
		double path = beam.ComputeIncidentOpticalPath(beam.direction, p);
		beams.internal.AddOpticalPath(path);
		beams.external.AddOpticalPath(path);
	}
	else
	{
		double path = beam.ComputeSegmentOpticalPath(reRiEff,
													 beams.internal.Center());
#ifdef _DEBUG // DEB
		beams.internal.ops = beam.ops;
		beams.external.ops = beam.ops;
#endif
		path += beam.opticalPath;
		beams.internal.AddOpticalPath(path);
		beams.external.AddOpticalPath(path);
	}
}

