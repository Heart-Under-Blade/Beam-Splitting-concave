#include <QtTest>

#include "ScatteringConvex.h"
#include "BulletRosette.h"
#include "Tracks.h"
#include "handler/HandlerPO.h"
#include "handler/HandlerPOTotal.h"
#include "global.h"

using namespace std;

class PO : public QObject
{
	Q_OBJECT

public:
	PO();
	~PO();

private slots:
	void test_Absorption();
	void test_Nan();

private:
	HandlerPO *h;
	Particle *pt;
	Scattering *sc;
	Light incidentLight;
};

PO::PO()
{
	incidentLight.direction = Point3f(0, 0, -1);
	incidentLight.polarizationBasis = Point3f(0, 1, 0);
	Point3f point = incidentLight.direction * pt->MaximalDimention()/2;
	incidentLight.direction.d_param = DotProduct(point, incidentLight.direction);

//	pt = new BulletRosette(complex(1.3116, 0), 42.04, 100,
//						   (42.04*sqrt(3)*tan(DegToRad(62)))/4);
	pt = new Particle();
	pt->SetFromFile("out15_mbs.dat");
	sc = new ScatteringConvex(pt, &incidentLight, true, 5);
	h = new HandlerPO(pt, &incidentLight, 1.6);
	ScatteringSphere sp(1, 64, 181);
	h->SetScattering(sc);
	h->SetScatteringSphere(sp);
	h->SetAbsorptionAccounting(true);
}

PO::~PO()
{
	delete pt;
	delete sc;
	delete h;
}

void PO::test_Absorption()
{
	pt->SetRefractiveIndex(complex(1.2893, 0));
	sc->SetSplitting(pt);

	vector<Beam> outBeams;
	pt->Rotate(179.34, 37, 0);
	sc->ScatterLight(179.34, 37, outBeams);

	QVERIFY(outBeams.size() == 434);

	for (unsigned i = 0; i < outBeams.size(); ++i)
	{
		Beam &beam = outBeams[i];

		if (beam.nActs > 0)
		{
			vector<int> tr;
			Tracks::RecoverTrack(beam, pt->nFacets, tr);

			double path = sc->MeasureFullOpticalPath(beam, beam.Center(), tr);
#ifdef _DEBUG // DEB
			if (fabs(path - beam.opticalPath) >= /*10e-4*/0.015)
				int ggg = 0;
#endif
			QVERIFY(fabs(path - beam.opticalPath) < /*10e-4*/0.015);
		}
	}
}

void PO::test_Nan()
{
	pt->SetRefractiveIndex(complex(1.2893, 3.5365e-4));
	sc->SetSplitting(pt);

	vector<Beam> outBeams;
	pt->Rotate(M_PI-6*0.15707963267948966,
			   M_PI+65*0.15707963267948966, 0);
	sc->ScatterLight(M_PI-6*0.15707963267948966,
					 M_PI+65*0.15707963267948966, outBeams);

	QVERIFY(outBeams.size() == 359);
	QVERIFY(outBeams[3].nVertices == 6);
	QVERIFY(outBeams[3].lastFacetId == 3);
	QVERIFY(outBeams[3].nActs == 0);
	QVERIFY(outBeams[3].front > -1.42 && outBeams[3].front < -1.41);

	for (Beam &beam : outBeams)
	{
		beam.polarizationBasis = beam.RotateSpherical(
					-incidentLight.direction, incidentLight.polarizationBasis);

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

QTEST_APPLESS_MAIN(PO)

#include "tst_po.moc"

//#include <QtTest/QtTest>

//class TestQString: public QObject
//{
//	Q_OBJECT
//private slots:
//	void toUpper();
//};

//void TestQString::toUpper()
//{
//	QString str = "Hello";
//	QCOMPARE(str.toUpper(), QString("HELLO"));
//}

//QTEST_APPLESS_MAIN(TestQString)

//#include "tst_po.moc"
