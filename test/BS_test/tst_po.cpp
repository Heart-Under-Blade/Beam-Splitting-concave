#include <QtTest>

#include "ScatteringNonConvex.h"
#include "BulletRosette.h"
#include "Tracks.h"
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

private:
	BulletRosette *pt;
	ScatteringNonConvex *sc;
	Light incidentLight;
};

PO::PO()
{
	incidentLight.direction = Point3f(0, 0, -1);
	incidentLight.polarizationBasis = Point3f(0, 1, 0);

	pt = new BulletRosette(complex(1.3116, 0), 42.04, 100,
						   (42.04*sqrt(3)*tan(DegToRad(62)))/4);
	sc = new ScatteringNonConvex(pt, &incidentLight, true, 4);
}

PO::~PO()
{
}

void PO::test_Absorption()
{
	Point3f point = incidentLight.direction * pt->GetRotationRadius();
	incidentLight.direction.d_param = DotProduct(point, incidentLight.direction);

	vector<Beam> outBeams;
	sc->ScatterLight(179.34, 37, outBeams);

	pt->Output();

//	QVERIFY(outBeams.size() == 3359);
//	Beam &b = outBeams[174];
//	QVERIFY(b.nActs == 1);
//	QVERIFY(b.size == 3);
//	QVERIFY(b.lastFacetId == 32);
//	QVERIFY(b.opticalPath >= 20022.8019);
//	QVERIFY(b.opticalPath <= 20022.8021);

	// absorption
	for (unsigned i = 0; i < outBeams.size(); ++i)
	{
		Beam &beam = outBeams[i];

		if (beam.nActs > 0)
		{
			vector<int> tr;
			Tracks::RecoverTrack(beam, pt->nFacets, tr);

			double path = sc->ComputeInternalOpticalPath(beam, tr);
#ifdef _DEBUG // DEB
			if (fabs(path - beam.opticalPath) >= 10e-4)
				int ggg = 0;
#endif
			QVERIFY(fabs(path - beam.opticalPath) < 10e-4);
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
