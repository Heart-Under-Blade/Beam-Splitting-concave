#include "Particle.h"
#include "common.h"
#include "geometry_lib.h"

#include <fstream>
#include <iostream>
#include <cstring>

Particle::Particle()
{
	SetFacetIndices();
}

Particle::Particle(int nFacets, bool isNonConvex)
	: m_isNonConvex(isNonConvex)
{
	nElems = nFacets;
	SetFacetIndices();
}

void Particle::SetFacetIndices()
{
	for (int i = 0; i < MAX_FACET_NUM; ++i)
	{
		elems[i].original.index = i;
		elems[i].actual.index = i;
	}
}

Facet *Particle::GetActualFacet(int i)
{
	return &elems[i].actual;
}

void Particle::SetFromFile(const std::string &filename, double reduceSize)
{
	std::ifstream pfile(filename, std::ios::in);

	if (!pfile.is_open())
	{
		std::cerr << "File \"" << filename << "\" is not found" << std::endl;
		throw std::exception();
	}

	const int bufSize = 1024;
	char *buff = (char*)malloc(sizeof(char) * bufSize);

	nElems = 0;
	Facet *facet = &(elems[nElems++].original);

	char *ptr, *trash;

	pfile.getline(buff, bufSize);
	m_isNonConvex = strtol(buff, &trash, 10);

	pfile.getline(buff, bufSize);
	isAggregated = strtol(buff, &trash, 10);

	// read symmetry params
	{
		pfile.getline(buff, bufSize);

		ptr = strtok(buff, " ");
		m_symmetry.zenith = Orientation::DegToRad(strtod(ptr, &trash));

		ptr = strtok(NULL, " ");
		m_symmetry.azimuth = Orientation::DegToRad(strtod(ptr, &trash));
	}

	pfile.getline(buff, bufSize); // skip empty line

	while (!pfile.eof()) // read vertices of facets
	{
#ifdef _DEBUG // DEB
		if (nElems == 35)
			int ff = 0;
#endif
		pfile.getline(buff, bufSize);
		ptr = strtok(buff, " ");

		if (strlen(buff) == 0)
		{
			facet = &(elems[nElems++].original);
			continue;
		}

		int c_i = 0;

		while (ptr != NULL)
		{
			double value = strtod(ptr, &trash);
			facet->vertices[facet->nVertices].coordinates[c_i++] = value;
			ptr = strtok(NULL, " ");
		}

		++(facet->nVertices);
	}

	pfile.close();

	// correction of number of facet
	while (elems[nElems-1].original.nVertices == 0)
	{
		--nElems;
	}

	if (reduceSize > 0)
	{
		ReduceSmallEdges(reduceSize);
	}

	SetDefaultNormals();
	ResetPosition();
	SetDefaultCenters();

	if (IsNonConvex() || isAggregated)
	{
		for (int i = 0; i < nElems; ++i)
		{
			elems[i].original.isOverlayedIn = true;
			elems[i].original.isOverlayedOut = true;
			elems[i].actual.isOverlayedIn = true;
			elems[i].actual.isOverlayedOut = true;
		}
	}
}

void Particle::Rotate(const Orientation &orientation)
{
	rotAngle = orientation;
	m_rotator.SetRotationAngle(rotAngle);

	// REF: слить всё в один цикл
	for (int i = 0; i < nElems; ++i)
	{
		auto &original = elems[i].original;
		auto &actual = elems[i].actual;

		for (int j = 0; j < original.nVertices; ++j)
		{
			m_rotator.RotatePoint(original.vertices[j], actual.vertices[j]);
		}

		// centers
		m_rotator.RotatePoint(original.center, actual.center);
	}

	RotateNormals();
}

void Particle::CommitState()
{
	for (int i = 0; i < nElems; ++i)
	{
		for (int j = 0; j < elems[i].actual.nVertices; ++j)
		{
			elems[i].original.vertices[j] = elems[i].actual.vertices[j];
		}
	}
}

void Particle::ProjectToXYPlane(std::vector<Point2f> &projection) const
{
	Point3f n(0, 0, 1);
	n.d_param = MaximalDimension();

	for (int i = 0; i < nElems; i++)
	{
		auto &f = elems[i].actual;

		if (Vector3f::IsAcuteAngled(f.in_normal, -n))
		{
			for (int j = 0; j < f.nVertices; j++)
			{
				double tmp = (n.d_param - Point3f::DotProduct(n, f.vertices[j]));
				auto p = f.vertices[j] + n*tmp;
				projection.push_back(Point2f(p.coordinates[0], p.coordinates[1]));
			}
		}
	}
}

void Particle::Concatenate(const std::vector<Particle> &parts)
{
	int i = 0;
	nElems = 0;

	for (const Particle &part : parts)
	{
		nElems += part.nElems;

		for (int j = 0; j < part.nElems; ++j)
		{
			elems[i++].original = part.elems[j].actual;
		}
	}

	isAggregated = true;
}

void Particle::RemoveFacet(int index)
{
	for (int i = index+1; i < nElems; ++i)
	{
		elems[i-1] = elems[i];
	}

	--nElems;

	SetFacetIndices();
}

double Particle::LongRadius() const
{
	Point3f p0(0, 0, 0);

	double radius = 0;

	for (int i = 0; i < nElems; ++i)
	{
		const Facet &facet = elems[i].actual;

		for (int j = 0; j < facet.nVertices; ++j)
		{
			Point3f v_len = facet.vertices[j] - p0;
			double len = Point3f::Length(v_len);

			if (len > radius)
			{
				radius = len;
			}
		}
	}

	return radius;
}

double Particle::Area() const
{
	double area = 0;

	for (int i = 0; i < nElems; ++i)
	{
		area += elems[i].original.Area();
	}

	return area;
}

Point3f Particle::Center() const
{
	Point3f center = Point3f(0, 0, 0);
	int nVertices = 0;

	for (int i = 0; i < nElems; ++i)
	{
		auto &facet = elems[i].original;
		nVertices += facet.nVertices;

		for (int j = 0; j < facet.nVertices; ++j)
		{
			center += facet.vertices[j];
		}
	}

	center /= nVertices;
	return center;
}

double Particle::Volume() const
{
	double volume = 0;
	Point3f center = Center();

	for (int i = 0; i < nElems; ++i)
	{
		const Facet &facet = elems[i].original;

		Point3f p = Geometry::ProjectPointToPlane(center, facet.ex_normal,
												  facet.in_normal);
		double h = Point3f::Length(p - center);
		volume += (facet.Area()*h)/3;
	}

	return volume;
}

void Particle::Scale(double ratio)
{
	for (int i = 0; i < nElems; ++i)
	{
		for (int j = 0; j < elems[i].original.nVertices; ++j)
		{
			elems[i].original.vertices[j] *= ratio;
		}
	}

	SetDefaultNormals();
	SetDParams();
	ResetPosition();
	SetDefaultCenters();
}

void Particle::Resize(double size)
{
	double dMax = MaximalDimension();
	double ratio = size/dMax;
	Scale(ratio);
	ResetPosition();
}

double Particle::MaximalDimension() const
{
	double Dmax = 0;
	double newDmax;

	Polygon pol;

	for (int i = 0; i < nElems; ++i)
	{
		pol.Concat(elems[i].original);
	}

	for (int i = 0; i < pol.nVertices; ++i)
	{
		for (int j = 0; j < pol.nVertices; ++j)
		{
			newDmax = Point3f::Length(pol.vertices[j] - pol.vertices[i]);

			if (newDmax > Dmax)
			{
				Dmax = newDmax;
			}
		}
	}

	return Dmax;
}

const Orientation &Particle::GetSymmetry() const
{
	return m_symmetry;
}

void Particle::GetFacets(int begin, int end, Array<Facet*> &facets)
{
	for (int i = begin; i < end; ++i)
	{
		Facet *f = &(elems[i].actual);
		facets.Add(f);
	}
}

void Particle::GetPartByFacet(Facet *facet, Array<Facet*> &facets)
{
	UNUSED(facet);
	GetFacets(0, nElems, facets);
}

void Particle::Move(float dx, float dy, float dz)
{
	for (int i = 0; i < nElems; ++i)
	{
		auto &origin = elems[i].original;
		auto &actual = elems[i].actual;

		for (int j = 0; j < origin.nVertices; ++j)
		{
			actual.vertices[j] = origin.vertices[j] + Point3f(dx, dy, dz);
		}
	}
}

bool Particle::IsNonConvex() const
{
	return m_isNonConvex;
}

void Particle::Output()
{
	std::ofstream file("particle.dat", std::ios::out);

	file << (int)m_isNonConvex << std::endl;
	file << (int)isAggregated << std::endl;
	file << Orientation::RadToDeg(m_symmetry.zenith) << ' '
		 << Orientation::RadToDeg(m_symmetry.azimuth) << std::endl << std::endl;

	for (int i = 0; i < nElems; ++i)
	{
		Facet &facet = elems[i].actual;

		for (int j = 0; j < facet.nVertices; ++j)
		{
			Point3f p = facet.vertices[j];
			file << p.coordinates[0]
					<< ' ' << p.coordinates[1]
					<< ' ' << p.coordinates[2];
			file << std::endl;
		}

		file << std::endl;
	}

	file.close();
}

int Particle::ReduceEdge(int facetNo, int i1, int i2)
{
	auto &p1 = elems[facetNo].original.vertices[i1];
	auto &p2 = elems[facetNo].original.vertices[i2];

	double eps = 0.001;
	Point3f newVertex = (p2 + p1)/2;

	int removedVertices = 0;
	int sameEdgeFacetNo = -1;

	for (int i = 0; i < nElems; ++i)
	{
		if (i != facetNo)
		{
			auto &facet = elems[i].original;
			int nFoundVertices = 0;

			for (int j = 0; j < facet.nVertices && nFoundVertices < 2; ++j)
			{
				auto &v = facet.vertices[j];

				if (v.IsEqualTo(p1, eps) || v.IsEqualTo(p2, eps))
				{
					++nFoundVertices;

					if (nFoundVertices < 2)
					{
						v = newVertex;
					}
					else
					{
						facet.RemoveVertex(j);
						++removedVertices;
						sameEdgeFacetNo = i;
					}
				}
			}
		}
	}

#ifdef _DEBUG // DEB
	if (removedVertices != 1)
	{
		int fff = 0;
	}
#endif

	elems[facetNo].original.vertices[i1] = newVertex;
	elems[facetNo].original.RemoveVertex(i2);

	return sameEdgeFacetNo;
}

void Particle::ReduceSmallEdges(double minSize)
{
	bool isReduced;

	do
	{
		isReduced = false;

		for (int i = 0; i < nElems; ++i)
		{
			auto &facet = elems[i].original;
			int jPrev = facet.nVertices-1;

			for (int j = 0; j < facet.nVertices; ++j)
			{
				Vector3f v = facet.vertices[j] - facet.vertices[jPrev];

				if (Point3f::Length(v) < minSize)
				{
#ifdef _DEBUG // DEB
					int nv = elems[i].original.nVertices;

					if (nv == 1)
						int ff = 1;
#endif
					int sameEdgeFacetNo = ReduceEdge(i, jPrev, j);
					isReduced = true;

					if (elems[i].original.nVertices < 3)
					{
						RemoveFacet(i);
					}

					if (sameEdgeFacetNo > 0)
					{
						if (elems[sameEdgeFacetNo].original.nVertices < 3)
						{
							RemoveFacet(sameEdgeFacetNo);
						}
					}

					break;
				}

				jPrev = j;
			}
		}
	} while (isReduced);
}

void Particle::SetDefaultNormals()
{
	for (int i = 0; i < nElems; ++i)
	{
		elems[i].original.SetNormal();
	}
}

void Particle::ResetPosition()
{
	for (int i = 0; i < nElems; ++i)
	{
		elems[i].actual = elems[i].original;
	}
}

void Particle::SetDefaultCenters()
{
	for (int i = 0; i < nElems; ++i)
	{
		elems[i].original.SetCenter();
	}
}

void Particle::RotateNormals()
{
	for (int i = 0; i < nElems; ++i)
	{
		auto &facet = elems[i];
		m_rotator.RotatePoint(facet.original.in_normal, facet.actual.in_normal);
	}

	SetDParams();

	for (int i = 0; i < nElems; ++i)
	{
		Facet &facet = elems[i].actual;
		facet.ex_normal = -facet.in_normal;
		facet.ex_normal.d_param = -facet.in_normal.d_param;
	}
}

void Particle::SetDParams()
{
	for (int i = 0; i < nElems; ++i)
	{
		Facet &facet = elems[i].actual;
		double d = Point3f::DotProduct(facet.vertices[0], facet.in_normal);
		facet.in_normal.d_param = -d;
	}
}

void Particle::SetSymmetry(double beta, double gamma)
{
	m_symmetry.zenith = beta;
	m_symmetry.azimuth = gamma;
}
