#pragma once

#include "types.h"
#include "PhysMtr.hpp"

class Tracer
{
public:
	Tracer();

private:
	Point3f m_startIncidentDir;
	Arr2D m_mxd;
};