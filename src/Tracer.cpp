#include "Tracer.h"

Tracer::Tracer()
	: m_mxd(0, 0, 0, 0)
{
	m_startIncidentDir = Point3f(0, 0, 1);
}
