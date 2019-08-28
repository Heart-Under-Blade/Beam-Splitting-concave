#pragma once

#include "Polygon.h"

class Facet : public Polygon
{
public:
	Point3f normal[2];	///< internal and external normals
	Point3f center;		///< center of facet polygon (for fast access without calc)

	bool isVisibleIn = true;
	bool isVisibleOut = true;

	void SetNormal();
	void SetCenter();

	Facet & operator = (const Facet &other);
};
