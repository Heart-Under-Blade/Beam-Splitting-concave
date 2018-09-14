#include "Incidence.h"

#include "Beam.h"
#include "Splitting.h"

void Incidence::ComputeOpticalPaths(const Beam &incidentBeam, Splitting &splitter) const
{
	if (incidentBeam.opticalPath < FLT_EPSILON)
	{
		Point3f p = splitter.inBeam.Center();
		double path = splitter.ComputeIncidentOpticalPath(incidentBeam.direction, p);
		splitter.inBeam.opticalPath = 0;
		splitter.outBeam.opticalPath = 0;
		splitter.inBeam.AddOpticalPath(path);
		splitter.outBeam.AddOpticalPath(path);
	}
	else
	{
		double path = splitter.ComputeSegmentOpticalPath(incidentBeam,
														 splitter.inBeam.Center());
		path += incidentBeam.opticalPath;
	#ifdef _DEBUG // DEB
		inBeam.ops = incidentBeam.ops;
		outBeam.ops = incidentBeam.ops;
		inBeam.ops.push_back(path);
		outBeam.ops.push_back(path);
	#endif
		splitter.inBeam.AddOpticalPath(path);
		splitter.outBeam.AddOpticalPath(path);
	}
}

