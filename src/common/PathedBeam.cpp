#include "PathedBeam.h"

PathedBeam::PathedBeam()
{
	opticalPath = 0;
}

void PathedBeam::Clear()
{
	Beam::Clear();
	opticalPath = 0;
}

void PathedBeam::AddOpticalPath(double path)
{
	opticalPath += path;
	front = Point3f::DotProduct(-direction, Center());
#ifdef _DEBUG // DEB
	ops.push_back(path);
#endif
}

double PathedBeam::ComputeOutgoingOpticalPath() const
{
	return FAR_ZONE_DISTANCE + front;
}

double PathedBeam::ComputeIncidentOpticalPath(const Point3f &direction,
											  const Point3f &facetPoint) const
{
	return FAR_ZONE_DISTANCE + Point3f::DotProduct(direction, facetPoint);
}

double PathedBeam::ComputeSegmentOpticalPath(const double &reRiEff,
											 const Point3f &facetPoint) const
{
	double tmp = Point3d::DotProduct(direction, facetPoint);
	double path = fabs(tmp + front); // refractive index of external media = 1

#ifdef _DEBUG // DEB
//	if (tmp + front < 0)
//		int fff = 0;
	Point3f dd = beam.Center();
	Point3f pd = beam.Center() + (beam.direction * 10);

	/* ПРоверка на нахождения плоскости за пучком
	 * REF: вынести в случай невыпуклых частиц, т.к. характерно только для них */
	double cosB = Point3d::DotProduct(Point3d::Normalize(facetPoint - beam.Center()), beam.direction);
#else
	Vector3f vectorToCenter = facetPoint - Center();
	Point3f::Normalize(vectorToCenter);
	double cosB = Point3f::DotProduct(vectorToCenter, direction);
#endif

	if (cosB < 0 && isInside)
	{
		path = -path;
	}

	if (isInside)
	{
		path *= sqrt(reRiEff);
	}

	return path;
}

PathedBeam &PathedBeam::operator =(const PathedBeam &other)
{
	if (this != &other) // OPT: попробовать убрать это уловие для ускорения
	{
		Beam::operator=(other);
		opticalPath = other.opticalPath;
		front = other.front;
	}

	return *this;
}

PathedBeam &PathedBeam::operator =(PathedBeam &&other)
{
	if (this != &other)
	{
		Beam::operator =(other);
		opticalPath = other.opticalPath;
		front = other.front;
		other.opticalPath = 0;
		other.front = 0;
#ifdef _DEBUG // DEB
		ops = other.ops;
#endif
	}

	return *this;
}
