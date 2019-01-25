#include "Converter.h"

#include <iomanip>
#include <float.h>
#include <fstream>
#include <cstring>
#include <limits.h>

#include "Particle.h"

#define NRM_EPS 0.01
#define PNT_EPS 10*FLT_EPSILON

bool IsNear(int i1, int i2, int n)
{
	return (abs(i2 - i1) == 1) ||
			(i2 == n-1 && i1 == 0) ||
			(i2 == 0 && i1 == n-1);
}

int FindEqualPoint(const Point3f &p, const Facet &merged)
{
	if (p.IsEqualTo(merged.arr[0], PNT_EPS))
	{
		return 0;
	}
	else if (p.IsEqualTo(merged.arr[1], PNT_EPS))
	{
		return 1;
	}
	else if (p.IsEqualTo(merged.arr[2], PNT_EPS))
	{
		return 2;
	}
	else // no equal points
	{
		return -1;
	}
}

bool FindEqualPoints(const Facet &origin, const Facet &next, Array<int> &points)
{
	int nextIndex;

	for (int i = 0; i < origin.nVertices && points.nElems != 6; ++i)
	{
		nextIndex = FindEqualPoint(origin.arr[i], next);

		if (nextIndex != -1)
		{
			points.Add(i);
			points.Add(nextIndex);
		}
	}

	return (points.nElems >= 4) &&
			IsNear(points.elems[0], points.elems[2], origin.nVertices);
}

Point3f Converter::ReadVertex(char *buff, char *ptr, char *trash)
{
	Point3f p;
	ptr = strtok(buff, " "); // skip "vertex" word

	ptr = strtok(NULL, " ");
	p.cx = strtod(ptr, &trash);
	ptr = strtok(NULL, " ");
	p.cy = strtod(ptr, &trash);
	ptr = strtok(NULL, " ");
	p.cz = strtod(ptr, &trash);

	return p;
}

Converter::Converter()
{

}

void Converter::ReadStl(const std::string &filename, std::vector<Facet> &triangles)
{
	std::ifstream pfile(filename, std::ios::in);

	if (!pfile.is_open())
	{
		std::cerr << "File \"" << filename << "\" is not found" << std::endl;
		throw std::exception();
	}

	const int bufSize = 1024;
	char *buff = (char*)malloc(sizeof(char) * bufSize);

	char *ptr, *trash;

	pfile.getline(buff, bufSize); // skip first line "solid HOLDER"
	pfile.getline(buff, bufSize);

	while (std::strncmp(buff, "endsolid", bufSize) < 0)
	{
		if (std::string(buff) == "")
		{
			pfile.getline(buff, bufSize); // skip empty line
			continue;
		}

		Facet facet;

		// read normal
		ptr = strtok(buff, " "); // skip "facet" word
		ptr = strtok(NULL, " "); // skip "normal" word

		ptr = strtok(NULL, " ");
		facet.ex_normal.cx = strtod(ptr, &trash);
		ptr = strtok(NULL, " ");
		facet.ex_normal.cy = strtod(ptr, &trash);
		ptr = strtok(NULL, " ");
		facet.ex_normal.cz = strtod(ptr, &trash);

		// read vertices
		pfile.getline(buff, bufSize); // skip "outer loop" line

		pfile.getline(buff, bufSize);
		facet.AddVertex(ReadVertex(buff, ptr, trash));
		pfile.getline(buff, bufSize);
		facet.AddVertex(ReadVertex(buff, ptr, trash));
		pfile.getline(buff, bufSize);
		facet.AddVertex(ReadVertex(buff, ptr, trash));

		triangles.push_back(facet);

		pfile.getline(buff, bufSize); // skip "endloop" line
		pfile.getline(buff, bufSize); // skip "endfacet" line
		pfile.getline(buff, bufSize);
	}

	pfile.close();
}

void Converter::ReadCry(const std::string &filename, std::vector<Facet> &crystal)
{
//    std::ifstream pfile(filename, std::ios::in);

//    if (!pfile.is_open())
//    {
//        std::cerr << "File \"" << filename << "\" is not found" << std::endl;
//        throw std::exception();
//    }

//    const int bufSize = 1024;
//    char *buff = (char*)malloc(sizeof(char) * bufSize);

//    char *ptr, *trash;

//    pfile.getline(buff, bufSize); // skip first line "solid HOLDER"
//    pfile.getline(buff, bufSize);

//    while (std::strncmp(buff, "endsolid", bufSize) < 0)
//    {
//        if (std::string(buff) == "")
//        {
//            pfile.getline(buff, bufSize); // skip empty line
//            continue;
//        }

//        Facet facet;

//        // read normal
//        ptr = strtok(buff, " "); // skip "facet" word
//        ptr = strtok(NULL, " "); // skip "normal" word

//        ptr = strtok(NULL, " ");
//        facet.ex_normal.cx = strtod(ptr, &trash);
//        ptr = strtok(NULL, " ");
//        facet.ex_normal.cy = strtod(ptr, &trash);
//        ptr = strtok(NULL, " ");
//        facet.ex_normal.cz = strtod(ptr, &trash);

//        // read vertices
//        pfile.getline(buff, bufSize); // skip "outer loop" line

//        pfile.getline(buff, bufSize);
//        facet.AddVertex(ReadVertex(buff, ptr, trash));
//        pfile.getline(buff, bufSize);
//        facet.AddVertex(ReadVertex(buff, ptr, trash));
//        pfile.getline(buff, bufSize);
//        facet.AddVertex(ReadVertex(buff, ptr, trash));

//        triangles.push_back(facet);

//        pfile.getline(buff, bufSize); // skip "endloop" line
//        pfile.getline(buff, bufSize); // skip "endfacet" line
//        pfile.getline(buff, bufSize);
//    }

//    pfile.close();
}

void Converter::MergeCrystal(std::vector<Facet> triangles,
							 std::vector<Facet> &mergedFacets)
{
	while (!triangles.empty())
	{
		std::vector<Facet> oneFacetTriangles;
		Point3f normal = triangles[0].Normal();
		double d1 = -Point3f::DotProduct(triangles[0].arr[0], -normal);

		for (int i = 0; i < triangles.size(); ++i)
		{
			const Facet &checking = triangles[i];

			if (normal.IsEqualTo(checking.Normal(), NRM_EPS))
			{
				double d2 = -Point3f::DotProduct(checking.arr[0], -normal);

				if (fabs(d1 - d2) < NRM_EPS)
				{
					oneFacetTriangles.push_back(checking);

					auto it = triangles.begin();
					std::advance(it, i);
					triangles.erase(it);
					--i;
				}
			}
		}

		while (!oneFacetTriangles.empty())
		{
			Facet merged;
			// OutputFacets(oneFacetTriangles);
			MergeTriangles(oneFacetTriangles, merged);
			mergedFacets.push_back(merged);
		}
	}
}

void Merge(const Array<int> &points, const Facet &checking, Facet &merged)
{
	if (points.nElems == 4)
	{
		int pointToInsert;

		if ((points.elems[1] == 0 && points.elems[3] == 1) ||
				(points.elems[1] == 1 && points.elems[3] == 0))
		{
			pointToInsert = 2;
		}
		else if ((points.elems[1] == 1 && points.elems[3] == 2) ||
				 (points.elems[1] == 2 && points.elems[3] == 1))
		{
			pointToInsert = 0;
		}
		else //
		{
			pointToInsert = 1;
		}

		int place = (points.elems[0] == merged.nVertices-1) ? 0
															: points.elems[0] + 1;

		if ((points.elems[0] == 0 && points.elems[2] == merged.nVertices-1) ||
				(points.elems[0] == merged.nVertices-1 && points.elems[2] == 0))
		{
			place = merged.nVertices;
		}

		merged.InsertVertex(place, checking.arr[pointToInsert]);
#ifdef _DEBUG // DEB
		if (merged.nVertices > 200)
			int gg = 0;
#endif
	}
	else if (points.nElems == 6)
	{
		if (IsNear(points.elems[0], points.elems[2], merged.nVertices) &&
				IsNear(points.elems[2], points.elems[4], merged.nVertices))
		{
			merged.DeleteVertex(points.elems[2]);
		}
	}
	else
	{
		throw std::exception();
	}
}

bool IsConvex(const Facet &merged)
{
	double dp1, dp2;
	int i1, i2;
	Vector3f v1, v2;

	Vector3f n = merged.Normal();

	v1 = merged.arr[1] - merged.arr[0];
	v2 = merged.arr[2] - merged.arr[0];
	dp1 = Vector3f::DotProduct(Vector3f::CrossProduct(v1, v2), n);

	for (int i = 1; i < merged.nVertices; ++i)
	{
		if (i == merged.nVertices-2)
		{
			i1 = merged.nVertices - 1;
			i2 = 0;
		}
		else if (i == merged.nVertices-1)
		{
			i1 = 0;
			i2 = 1;
		}
		else
		{
			i1 = i + 1;
			i2 = i + 2;
		}

		v1 = merged.arr[i1] - merged.arr[i];
		v2 = merged.arr[i2] - merged.arr[i];
		dp2 = Vector3f::DotProduct(Vector3f::CrossProduct(v1, v2), n);

		if (dp2*dp1 < 0)
		{
			return false;
		}
	}

	return true;
}

void Converter::MergeTriangles(std::vector<Facet> &rest, Facet &convex)
{
	std::vector<Facet> triangles = rest;
	Facet merged = triangles[0];
#ifdef _DEBUG // DEB
	if (triangles.size() == 88)
		int fff = 0;
#endif
	auto it = triangles.begin();
	std::advance(it, 0);
	triangles.erase(it);

	convex = merged;
	rest = triangles;

	while (!triangles.empty())
	{
		bool isAdded = false;

		for (int i = 0; i < triangles.size(); ++i)
		{
			const Facet &triangle = triangles[i];
			Array<int> points;

			if (FindEqualPoints(merged, triangle, points))
			{
				Merge(points, triangle, merged);

				it = triangles.begin();
				std::advance(it, i);
				triangles.erase(it);

				if (IsConvex(merged))
				{
					convex = merged;
					rest = triangles;
				}

				isAdded = true;
			}
		}

		if (!isAdded)
		{
			break;
		}
	}
}

void Converter::WriteStl(const std::vector<Facet> &triangles,
						 const std::string &outFile)
{
	std::ofstream ofile(outFile + ".stl", std::ios::out);

	if (!ofile.is_open())
	{
		std::cerr << "File \"" << outFile << "\" is not found" << std::endl;
		throw std::exception();
	}

	ofile << std::setprecision(6) << std::scientific;

	int offset = 3;

	int nSpaces = offset;
	ofile << "solid HOLDER" << std::endl;

	for (const Facet &facet : triangles)
	{
		const Point3f &n = facet.Normal();
		ofile << std::string(nSpaces, ' ') << "facet normal "
			  << n.point[0] << ' '
			  << n.point[1] << ' '
			  << n.point[2] << std::endl;

		nSpaces += offset;
		ofile << std::string(nSpaces, ' ') << "outer loop" << std::endl;

		nSpaces += offset;

		for (int i = 0; i < facet.nVertices; ++i)
		{
			ofile << std::string(nSpaces, ' ') << "vertex "
				  << facet.arr[i].point[0] << ' '
				  << facet.arr[i].point[1] << ' '
				  << facet.arr[i].point[2] << std::endl;
		}

		nSpaces -= offset;
		ofile << std::string(nSpaces, ' ') << "endloop" << std::endl;

		nSpaces -= offset;
		ofile << std::string(nSpaces, ' ') << "endfacet" << std::endl;
	}

	ofile << "endsolid" << std::endl;
	ofile.close();
}

void Converter::WriteNat(const std::vector<Facet> &crystal,
						 const std::string &outFile)
{
	std::string crystalName = outFile;
	std::ofstream ofile(outFile + ".dat", std::ios::out);

	if (!ofile.is_open())
	{
		std::cerr << "File \"" << outFile << "\" is not found" << std::endl;
		throw std::exception();
	}

	ofile << std::setprecision(10);

	int maxNVertices = 0;
	int maxNVertInFacet = 0;

	for (const Facet &f : crystal)
	{
		maxNVertices += f.nVertices;

		if (f.nVertices > maxNVertInFacet)
		{
			maxNVertInFacet = f.nVertices;
		}
	}

	ofile << "class " << crystalName <<  ": public Crystal {\n\
	double tt, ps, fi; // orientation of prizm\n\
	double w = 100;       // sizes of prizm: radius and semiheight\n\
	double __fastcall Second(void) { return this->w; }\n\
	double __fastcall Third(void) { throw \" Prizm::Third(): Error! \"; }\n\
	double __fastcall Forth(void) { throw \" Prizm::Forth(): Error! \"; }\n\
public:\n\
	// constructor and destructor\n\
	__fastcall " << crystalName <<  "(const complex& r) :\n\
	" << crystalName <<  "(r," << maxNVertices << ',' << crystal.size() << ','
						<< maxNVertInFacet << "), tt(0), ps(0), fi(0) \n\
	{ this->SetVertices(); this->SetFacets(); }\n\
	virtual __fastcall ~TypeCrystal(void) {}\n\
	// parameters\n\
	double __fastcall  First(void)  { return this->w; }\n\
	void   __fastcall  SetVertices(void); // sets vertices of prizm\n\
	void   __fastcall  SetFacets(void); // sets facets of prizm\n\
};\n\
\n\
\n\
void __fastcall " << crystalName <<  "::SetVertices(void)\n\
{\n";

	 int vertexCount = 0;

	 for (int i = 0; i < crystal.size(); ++i)
	 {
		 for (int j = 0; j < crystal[i].nVertices; ++j)
		 {
			 const Point3f &p = crystal[i].arr[j];
			 ofile << "\tthis->p[" << vertexCount++ << "] = Point3d(" << p.point[0] << ", " << p.point[1] << ", " << p.point[2] << ");" << std::endl;
		 }
	 }

	 ofile << "Point3d cnt = CenterOfGravity();\n\
			  for (int i=0; i<this->M; i++)\n\
			  this->p[i] -= cnt;\n\
}\n\
\n\
\n\
void __fastcall " << crystalName <<  "::SetFacets(void)\n\
{\n\
\tint* q;   // the basal facets:";

	vertexCount = 0;

	for (int i = 0; i < crystal.size(); ++i)
	{
		ofile << "\n\tq=this->Gr[" << i << "]; ";
		int first = vertexCount;

		for (int j = 0; j <= maxNVertInFacet; ++j)
		{
			if (j == crystal[i].nVertices)
			{
				ofile << "q[" << j << "]=" << first << "; ";
			}
			else if (j > crystal[i].nVertices)
			{
				ofile << "q[" << j << "]=" << -1 << "; ";
			}
			else
			{
				ofile << "q[" << j << "]=" << vertexCount++ << "; ";
			}
		}
	}

	ofile << "\n}\n";

	ofile.close();
}
