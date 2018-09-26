#include "Incidence.h"

#include "Beam.h"
#include "Splitting.h"

void Incidence::ComputeOpticalPaths(const Beam &beam, Splitting &splitter) const
{
	if (beam.opticalPath < FLT_EPSILON)
	{
		Point3f p = splitter.inBeam.Center();
		double path = splitter.ComputeIncidentOpticalPath(beam.direction, p);
		splitter.inBeam.opticalPath = 0;
		splitter.outBeam.opticalPath = 0;
		splitter.inBeam.AddOpticalPath(path);
		splitter.outBeam.AddOpticalPath(path);
	}
	else
	{
		double path = splitter.ComputeSegmentOpticalPath(beam, splitter.inBeam.Center());
#ifdef _DEBUG // DEB
		splitter.inBeam.ops = beam.ops;
		splitter.outBeam.ops = beam.ops;
		splitter.inBeam.ops.push_back(path);
		splitter.outBeam.ops.push_back(path);
#endif
		path += beam.opticalPath;
		splitter.inBeam.AddOpticalPath(path);
		splitter.outBeam.AddOpticalPath(path);
	}
}

