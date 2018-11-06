#include "Incidence.h"

void Incidence::ComputeOpticalPaths(const PathedBeam &beam, SplittedBeams<PathedBeam> &beams,
									Splitting &splitter) const
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
		double path = beam.ComputeSegmentOpticalPath(splitter.reRiEff,
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

