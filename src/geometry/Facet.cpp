#include "Facet.h"

Facet::Facet()
{
}

void Facet::SetNormal()
{
	ex_normal = Normal();
	in_normal = -ex_normal;
}

void Facet::SetCenter()
{
	center = Center();
}

Facet &Facet::operator =(const Facet &other)
{
	if (this != &other)
	{
		Polygon1::operator =(other);
		in_normal = other.in_normal;
		ex_normal = other.ex_normal;
		center = other.center;
	}

	return *this;
}

