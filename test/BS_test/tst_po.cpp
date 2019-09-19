#include <QtTest>

#include "ScatteringConvex.h"
#include "RegularColumn.h"
#include "Tracks.h"
#include "handler/HandlerPO.h"

using namespace std;

class PO : public QObject
{
	Q_OBJECT

public:
	PO();
	~PO();

private slots:
	void test_Absorption();
	void test_Particle();
	void test_BenchmarkScatter();
	void test_Nan();
	void test_NanCol();

private:
	HandlerPO *h;
	Particle *pt;
	Scattering *sc;
	Tracks *trs;
};

PO::PO()
{
//	pt = new BulletRosette(complex(1.3116, 0), 42.04, 100,
//						   (42.04*sqrt(3)*tan(DegToRad(62)))/4);
	pt = new Particle();
	pt->SetFromFile("out15_mbs.dat");

	sc = new ScatteringConvex(complex(), 5);
	sc->SetActiveParticle(pt);

	h = new HandlerPO(sc, 1.6);
	ScatteringSphere sp(1, 64, 181);
	h->SetScatteringSphere(sp);
	h->EnableAbsorption(true);
}

PO::~PO()
{
	delete pt;
	delete sc;
	delete h;
	delete trs;
}

void PO::test_Absorption()
{
	sc->SetRefractiveIndex(complex(1.2893, 0));

	vector<Beam> outBeams;

	pt->Rotate(Orientation(179.34, 37));
	sc->Scatter(&outBeams);

	QVERIFY(outBeams.size() == 434);

	// absorption
//	for (unsigned i = 0; i < outBeams.size(); ++i)
//	{
//		Beam &beam = outBeams[i];
//		BeamInfo info = h->ComputeBeamInfo(beam);

//		if (beam.nActs > 0)
//		{
//			vector<int> tr;
//			h->m_tracks->RecoverTrack(beam, tr);

//			double path = sc->MeasureFullOpticalPath(info, beam.Center());
//			QVERIFY(fabs(path - beam.opticalPath) < /*10e-4*/0.015);
//		}
//	}
}

void PO::test_Particle()
{
	// test normals
	for (int i = 0; i < pt->nElems; ++i)
	{
		Facet *f = pt->GetActualFacet(i);
		Point3f v;
		QVERIFY(fabs(Point3f::Length(f->normal[0]) - 1) < FLT_EPSILON);
		QVERIFY(fabs(Point3f::Length(f->normal[1]) - 1) < FLT_EPSILON);
	}
}

void PO::test_BenchmarkScatter()
{
	vector<Beam> outBeams;

	QBENCHMARK
	{
		for (double z = 0; z < 179.34; z += 2.1)
		{
			pt->Rotate(Orientation(z, 37));
			sc->Scatter(&outBeams);
		}
	}
}

void PO::test_Nan()
{
	trs = new Tracks(pt->nElems);
	h->SetTracks(trs);

	sc->SetRefractiveIndex(complex(1.2893, 3.5365e-4));

	vector<Beam> outBeams;
	pt->Rotate(Orientation(M_PI-6*0.15707963267948966,
						   M_PI+65*0.15707963267948966));
	sc->Scatter(&outBeams);

	QVERIFY(outBeams.size() == 359);
	QVERIFY(outBeams[3].nVertices == 6);
	QVERIFY(outBeams[3].facet->index == 3);
	QVERIFY(outBeams[3].actNo == 0);
//	QVERIFY(outBeams[3].front > -1.42 && outBeams[3].front < -1.41);

	const Beam &startBeam = sc->GetStartBeam();

	for (Beam &beam : outBeams)
	{
		beam.polarizationBasis = beam.RotateSpherical(
					-startBeam.direction, startBeam.polarizationBasis);

		BeamInfo info = h->ComputeBeamInfo(beam);

		if (h->m_isBadBeam)
		{
			continue;
		}

		for (int i = 0; i <= h->m_sphere.nAzimuth; ++i)
		{
			for (int j = 0; j <= h->m_sphere.nZenith; ++j)
			{
				Point3d &dir = h->m_sphere.directions[i][j];
				Point3d &vf = (j == 0) ? h->m_sphere.vf.back() : h->m_sphere.vf[i];
				matrixC diffractedMatrix = h->ApplyDiffraction(beam, info, dir, vf);

				double ddd = real(diffractedMatrix[0][0]);
#ifdef _DEBUG // DEB
				if (isnan(ddd))
					int fff = 0;
#endif
				QVERIFY(!isnan(ddd));
			}
		}
	}
}

void PO::test_NanCol()
{
	delete pt;
	pt = new RegularColumn(Size(68.9, 100));
	sc->SetRefractiveIndex(complex(1.2893, 3.5365e-4));

	delete trs;
	trs = new Tracks(pt->nElems);
	h->SetTracks(trs);

	vector<Beam> outBeams;
	pt->Rotate(Orientation(M_PI-0.15707963267948966,
						   M_PI+1.5707963267948966));
	sc->Scatter(&outBeams);
	const Beam &startBeam = sc->GetStartBeam();
#ifdef _DEBUG // DEB
	int cc = 0;
#endif
	for (Beam &beam : outBeams)
	{
#ifdef _DEBUG // DEB
		++cc;
		if (cc == 5)
			int ffffd = 0;
#endif
		h->m_isBadBeam = false;

		beam.polarizationBasis = beam.RotateSpherical(
					-startBeam.direction,
					startBeam.polarizationBasis);

		BeamInfo info = h->ComputeBeamInfo(beam);

		if (h->m_isBadBeam)
		{
			continue;
		}

		for (int i = 0; i <= h->m_sphere.nAzimuth; ++i)
		{
			for (int j = 0; j <= h->m_sphere.nZenith; ++j)
			{
				Point3d &dir = h->m_sphere.directions[i][j];
				Point3d &vf = (j == 0) ? h->m_sphere.vf.back() : h->m_sphere.vf[i];
				matrixC diffractedMatrix = h->ApplyDiffraction(beam, info, dir, vf);

				complex fff = diffractedMatrix[0][0];
				if (isnan(real(fff)))
					int ffff = 0;
				QVERIFY(!isnan(real(fff)));
			}
		}
	}
}

QTEST_APPLESS_MAIN(PO)

#include "tst_po.moc"
