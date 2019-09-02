#include "Particle.h"
#include <fstream>
#include <iostream>
#include <cstring>
#include "global.h"

Particle::Particle()
{
<<<<<<< HEAD
	isConcave = false;
=======
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
>>>>>>> origin/refactor
}

void Particle::SetFromFile(const std::string &filename)
{
	std::ifstream pfile(filename, std::ios::in);

	if (!pfile.is_open())
	{
		std::cerr << "File \"" << filename << "\" is not found" << std::endl;
		throw std::exception();
	}

	const int bufSize = 1024;
	char *buff = (char*)malloc(sizeof(char) * bufSize);

	nFacets = 0;
	Facet *facet = &(defaultFacets[nFacets++]);

	char *ptr, *trash;

	pfile.getline(buff, bufSize);
	isConcave = strtol(buff, &trash, 10);

	pfile.getline(buff, bufSize);
	isAggregated = strtol(buff, &trash, 10);

	// read symmetry params
	{
		pfile.getline(buff, bufSize);

		ptr = strtok(buff, " ");
		m_symmetry.beta = DegToRad(strtod(ptr, &trash));

		ptr = strtok(NULL, " ");
		m_symmetry.gamma = DegToRad(strtod(ptr, &trash));
	}

	pfile.getline(buff, bufSize); // skip empty line

	while (!pfile.eof()) // read vertices of facets
	{
		pfile.getline(buff, bufSize);
		ptr = strtok(buff, " ");

		if (strlen(buff) == 0)
		{
			facet = &(defaultFacets[nFacets++]);
			continue;
		}

		int c_i = 0;

		while (ptr != NULL)
		{
<<<<<<< HEAD
			facet->arr[facet->nVertices].coordinates[c_i++] = strtod(ptr, &trash);
=======
			double value = strtod(ptr, &trash);
			facet->vertices[facet->nVertices].coordinates[c_i++] = value;
>>>>>>> origin/refactor
			ptr = strtok(NULL, " ");
		}

		++(facet->nVertices);
	}

	pfile.close();

	// correction of number of facet
<<<<<<< HEAD
	while (defaultFacets[nFacets-1].nVertices == 0)
=======
	while (elems[nElems-1].original.nVertices == 0)
	{
		--nElems;
	}

	if (reduceSize > 0)
>>>>>>> origin/refactor
	{
		--nFacets;
	}

	SetDefaultNormals();
	Reset();
	SetDefaultCenters();

	if (isConcave || isAggregated)
	{
		for (int i = 0; i < nFacets; ++i)
		{
			defaultFacets[i].isVisibleIn = false;
			defaultFacets[i].isVisibleOut = false;
			facets[i].isVisibleIn = false;
			facets[i].isVisibleOut = false;
		}
	}
}

void Particle::Init(int facetCount, const complex &refrIndex)
{
	nFacets = facetCount;
	m_refractiveIndex = refrIndex;
}

void Particle::RotateCenters()
{
	for (int i = 0; i < nFacets; ++i)
	{
		RotatePoint(defaultFacets[i].center, facets[i].center);
	}
}

void Particle::Rotate(double beta, double gamma, double alpha)
{
	rotAngle = Angle{alpha, beta, gamma};
	SetRotateMatrix(beta, gamma, alpha);

	// REF: слить всё в один цикл
	for (int i = 0; i < nFacets; ++i)
	{
<<<<<<< HEAD
		for (int j = 0; j < facets[i].nVertices; ++j)
=======
		auto &original = elems[i].original;
		auto &actual = elems[i].actual;

		for (int j = 0; j < original.nVertices; ++j)
>>>>>>> origin/refactor
		{
			RotatePoint(defaultFacets[i].arr[j], facets[i].arr[j]);
		}
	}

	RotateNormals();
	RotateCenters();
}

<<<<<<< HEAD
void Particle::Fix()
{
	for (int i = 0; i < nFacets; ++i)
	{
		for (int j = 0; j < facets[i].nVertices; ++j)
		{
			defaultFacets[i].arr[j] = facets[i].arr[j];
		}
	}
}

void Particle::Scale(double ratio)
=======
void Particle::CommitState()
>>>>>>> origin/refactor
{
	for (int i = 0; i < nFacets; ++i)
	{
<<<<<<< HEAD
		for (int j = 0; j < defaultFacets[i].nVertices; ++j)
=======
		for (int j = 0; j < elems[i].actual.nVertices; ++j)
>>>>>>> origin/refactor
		{
			defaultFacets[i].arr[j].cx *= ratio;
			defaultFacets[i].arr[j].cy *= ratio;
			defaultFacets[i].arr[j].cz *= ratio;
		}
	}
}

<<<<<<< HEAD
void Particle::Resize(double size)
{
	double dMax = MaximalDimention();
	double ratio = size/dMax;
	Scale(ratio);
	Reset();
}

=======
>>>>>>> origin/refactor
void Particle::Concate(const std::vector<Particle> &parts)
{
	int i = 0;
	nFacets = 0;

	for (const Particle &part : parts)
	{
		nFacets += part.nFacets;

		for (int j = 0; j < part.nFacets; ++j)
		{
			defaultFacets[i++] = part.facets[j];
		}
	}

	isAggregated = true;
}

<<<<<<< HEAD
=======
void Particle::RemoveFacet(int index)
{
	for (int i = index+1; i < nElems; ++i)
	{
		elems[i-1] = elems[i];
	}

	--nElems;

	SetFacetIndices();
}

>>>>>>> origin/refactor
double Particle::LongRadius() const
{
	Point3f p0(0, 0, 0);

	double radius = 0;

	for (int i = 0; i < nFacets; ++i)
	{
<<<<<<< HEAD
		for (int j = 0; j < facets[i].nVertices; ++j)
=======
		const Facet &facet = elems[i].actual;

		for (int j = 0; j < facet.nVertices; ++j)
>>>>>>> origin/refactor
		{
			Point3f v_len = facets[i].arr[j] - p0;
			double len = Length(v_len);

			if (len > radius)
			{
				radius = len;
			}
		}
	}

	return radius;
}

<<<<<<< HEAD
double Particle::MaximalDimention() const
{
	double Dmax = 0;
	double newDmax;

	Polygon512 pol;

	for (int i = 0; i < nFacets; ++i)
	{
		pol.Concat(defaultFacets[i]);
	}

	for (int i = 0; i < pol.nVertices; ++i)
	{
		for (int j = 0; j < pol.nVertices; ++j)
		{
			newDmax = Length(pol.arr[j] - pol.arr[i]);

			if (newDmax > Dmax)
			{
				Dmax = newDmax;
			}
		}
	}

	return Dmax;
}

double Particle::Area()
=======
double Particle::Area() const
>>>>>>> origin/refactor
{
	double area = 0;

	for (int i = 0; i < nFacets; ++i)
	{
		area += facets[i].Area();
	}

	return area;
}

<<<<<<< HEAD
const complex &Particle::GetRefractiveIndex() const
=======
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
>>>>>>> origin/refactor
{
	return m_symmetry;
}

<<<<<<< HEAD
const Symmetry &Particle::GetSymmetry() const
=======
void Particle::GetFacets(int end, int begin, Array<Facet*> &facets)
>>>>>>> origin/refactor
{
	for (int i = begin; i < end; ++i)
	{
		Facet *f = &(elems[i].actual);
		facets.Add(f);
	}
}

void Particle::GetPartByFacet(Facet */*facet*/, Array<Facet*> &facets)
{
	GetFacets(0, nElems, facets);
}

void Particle::Move(float dx, float dy, float dz)
{
	for (int i = 0; i < nFacets; ++i)
	{
<<<<<<< HEAD
		for (int j = 0; j < defaultFacets[i].nVertices; ++j)
=======
		auto &origin = elems[i].original;
		auto &actual = elems[i].actual;

		for (int j = 0; j < origin.nVertices; ++j)
>>>>>>> origin/refactor
		{
			facets[i].arr[j] = defaultFacets[i].arr[j] + Point3f(dx, dy, dz);
		}
	}
}

bool Particle::IsConcave() const
{
	return isConcave;
}

void Particle::Output()
{
	std::ofstream M("particle.dat", std::ios::out);

	for (int i = 0; i < nFacets; ++i)
	{
<<<<<<< HEAD
		for (int j = 0; j < facets[i].nVertices; ++j)
		{
			Point3f p = facets[i].arr[j];
			M << p.coordinates[0] << ' '
							<< p.coordinates[1] << ' '
							<< p.coordinates[2] << ' '
							<< i ;
=======
		Facet &facet = elems[i].actual;

		for (int j = 0; j < facet.nVertices; ++j)
		{
			Point3f p = facet.vertices[j];
			M << p.coordinates[0]
					<< ' ' << p.coordinates[1]
					<< ' ' << p.coordinates[2]
					<< ' ' << i ;
>>>>>>> origin/refactor
			M << std::endl;
		}

		M << std::endl << std::endl;;
	}

	M.close();
}

<<<<<<< HEAD
void Particle::SetRefractiveIndex(const complex &value)
{
	m_refractiveIndex = value;
=======
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
>>>>>>> origin/refactor
}

void Particle::SetDefaultNormals()
{
	for (int i = 0; i < nFacets; ++i)
	{
		defaultFacets[i].SetNormal();
	}
}

void Particle::Reset()
{
	for (int i = 0; i < nFacets; ++i)
	{
		facets[i] = defaultFacets[i];
	}
}

void Particle::SetDefaultCenters()
{
	for (int i = 0; i < nFacets; ++i)
	{
		defaultFacets[i].SetCenter();
	}
}

void Particle::SetRotateMatrix(double beta, double gamma, double alpha)
{
	double cosA, cosB, cosG,
			sinA, sinB, sinG;

	sincos(alpha, &sinA, &cosA);
	sincos(beta,  &sinB, &cosB);
	sincos(gamma, &sinG, &cosG);

	double cosAcosB = cosA*cosB;
	double sinAcosG = sinA*cosG;
	double sinAsinG = sinA*sinG;

	m_rotMatrix[0][0] = cosAcosB*cosG - sinAsinG;
	m_rotMatrix[1][0] = sinAcosG*cosB + cosA*sinG;
	m_rotMatrix[2][0] = -sinB*cosG;

	m_rotMatrix[0][1] = -(cosAcosB*sinG + sinAcosG);
	m_rotMatrix[1][1] = cosA*cosG - sinAsinG*cosB;
	m_rotMatrix[2][1] = sinB*sinG;

	m_rotMatrix[0][2] = cosA*sinB;
	m_rotMatrix[1][2] = sinA*sinB;
	m_rotMatrix[2][2] = cosB;
}

void Particle::RotateNormals()
{
	for (int i = 0; i < nFacets; ++i)
	{
		RotatePoint(defaultFacets[i].in_normal, facets[i].in_normal);
	}

	SetDParams();

	for (int i = 0; i < nFacets; ++i)
	{
		facets[i].ex_normal = -facets[i].in_normal;
		facets[i].ex_normal.d_param = -facets[i].in_normal.d_param;
	}
}

void Particle::SetDParams()
{
	for (int i = 0; i < nFacets; ++i)
	{
		double d = DotProduct(facets[i].arr[0], facets[i].in_normal);
		facets[i].in_normal.d_param = -d;
	}
}

void Particle::RotatePoint(const Point3f &point, Point3f &result)
{
	result.cx = point.cx*m_rotMatrix[0][0] + point.cy*m_rotMatrix[0][1] + point.cz*m_rotMatrix[0][2];
	result.cy = point.cx*m_rotMatrix[1][0] + point.cy*m_rotMatrix[1][1] + point.cz*m_rotMatrix[1][2];
	result.cz = point.cx*m_rotMatrix[2][0] + point.cy*m_rotMatrix[2][1] + point.cz*m_rotMatrix[2][2];
}

void Particle::SetSymmetry(double beta, double gamma, double alpha)
{
	m_symmetry.beta = beta;
	m_symmetry.gamma = gamma;
	m_symmetry.alpha = alpha;
}
