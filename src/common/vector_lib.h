#pragma once

class Point3f;


double DotProduct(const Point3f &v1, const Point3f &v2);

double Norm(const Point3f &point);

void CrossProduct(const Point3f &v1, const Point3f &v2, Point3f &res);

void Normalize(Point3f &v);

Point3f NormalToFacet(const Point3f *facet);

void CopyPoints(Point3f *points, Point3f *result, int size);

Point3f CenterOfFacet(const Point3f *facet, int size);
