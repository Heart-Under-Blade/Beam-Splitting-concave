#include <QtTest>

#include "ScatteringNonConvex.h"
#include "BulletRosette.h"
#include "HollowColumn.h"
#include "Tracks.h"
#include "common.h"

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
	Particle *pt;
	ScatteringNonConvex *sc;
	Light incidentLight;
};

PO::PO()
{
	incidentLight.direction = Point3f(0, 0, -1);
	incidentLight.polarizationBasis = Point3f(0, 1, 0);
//	pt = new BulletRosette(complex(1.3116, 0), Size(42.04, 100),
//						   (42.04*sqrt(3)*tan(Angle::DegToRad(62)))/4);
	pt = new HollowColumn(complex(1.3116, 0), Size(42.04, 100), 5);
	sc = new ScatteringNonConvex(pt, incidentLight, 4);
}

PO::~PO()
{
}

void PO::test_Absorption()
{
	Point3f point = incidentLight.direction * pt->ComputeRotationRadius();
	incidentLight.direction.d_param = Point3f::DotProduct(point, incidentLight.direction);

	vector<Beam> outBeams;
	pt->Rotate(Angle3d(0, Angle3d::DegToRad(179.34), Angle3d::DegToRad(37)));
//	sc->RotateParticle(Angle(0, 179.34, 37));
	sc->ScatterLight(outBeams);

	pt->Output();

//	QVERIFY(outBeams.size() == 3359);
//	Beam &b = outBeams[174];
//	QVERIFY(b.nActs == 1);
//	QVERIFY(b.size == 3);
//	QVERIFY(b.lastFacetId == 32);
//	QVERIFY(b.opticalPath >= 20022.8019);
//	QVERIFY(b.opticalPath <= 20022.8021);

	QVERIFY(!outBeams.empty());

	// absorption

	Tracks tracks(pt->nElems);
	bool isCatchedBeams = false;

	for (unsigned i = 0; i < outBeams.size(); ++i)
	{
		Beam &beam = outBeams[i];
//		QVERIFY(!isnan(beam.front));

		if (beam.actNo > 0)
		{
			isCatchedBeams = true;
			vector<int> tr;
			tracks.RecoverTrack(beam, tr);

			OpticalPath path = sc->ComputeOpticalPath(beam, beam.Center(), tr);
			double total = path.GetTotal();
#ifdef _DEBUG // DEB
//			if (fabs(total - beam.opticalPath) >= 10e-4)
//				int ggg = 0;
//			QVERIFY(fabs(total - beam.opticalPath) < 10e-4);
#endif
		}
	}

	QVERIFY(isCatchedBeams);
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
