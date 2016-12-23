#pragma once

class Point3f;
class Point3d;
class Polygon;

float DotProduct(const Point3f &v1, const Point3f &v2);

double DotProductD(const Point3d &v1, const Point3d &v2);

double Norm(const Point3f &point);

void CrossProduct(const Point3f &v1, const Point3f &v2, Point3f &res);

void Normalize(Point3f &v);

Point3f NormalToPolygon(const Point3f *facet);

Point3f CenterOfPolygon(const Polygon &polygon);

double Length(const Point3f &v);
