TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

DESTDIR = ../bin

QMAKE_CXXFLAGS += -std=gnu++11
QMAKE_CXXFLAGS += -march=corei7 -msse4.2

CONFIG(release, debug|release): DEFINES += _NDEBUG
CONFIG(debug,	debug|release): DEFINES += _DEBUG

VERSION = 1.0.0

INCLUDEPATH += \
	../src \
	../src/bigint \
	../src/math \
	../src/particle \
	../src/geometry \
	../src/geometry/clipper \
	../src/common \
	../src/tracing

HEADERS += \
	../src/Beam.h \
	../src/geometry/Intersection.h \
	../src/test.h \
	../src/geometry/intrinsic/intrinsics.h \
    ../src/math/compl.hpp \
    ../src/math/matrix.hpp \
	../src/math/service.hpp \
    ../src/particle/Hexagonal.h \
    ../src/particle/Particle.h \
	../src/particle/ConcaveHexagonal.h \
    ../src/math/JonesMatrix.h \
    ../src/math/Mueller.hpp \
    ../src/math/PhysMtr.hpp \
    ../src/tracing/Tracing.h \
	../src/common/global.h \
	../src/common/macro.h \
    ../src/tracing/TracingConvex.h \
    ../src/tracing/TracingConcave.h \
	../src/CalcTimer.h \
	../src/Tracer.h \
	../src/geometry/clipper/clipper.hpp \
	../src/geometry/clipper/BeamClipper.h \
	../src/geometry/geometry_lib.h \
	../src/PhisBeam.h \
	../src/particle/TiltedHexagonal.h \
	../src/argparse.hpp \
    ../src/ArgParser.h \
    ../src/particle/HexagonalAggregate.h \
	../src/particle/CertainAggregate.h\
	../src/ArgPP.h \
    ../src/geometry/Polygon.h \
    ../src/geometry/Facet.h \
    ../src/particle/BulletRosette.h \
    ../src/particle/Bullet.h \
    ../src/common/Matrix4x4.h \
    ../src/common/MullerMatrix.h

SOURCES += \
    ../src/Beam.cpp \
	../src/geometry/Intersection.cpp \
	../src/main.cpp \
	../src/geometry/intrinsic/intrinsics.cpp \
    ../src/math/compl.cpp \
	../src/math/matrix.cpp \
    ../src/particle/Hexagonal.cpp \
	../src/particle/ConcaveHexagonal.cpp \
    ../src/particle/Particle.cpp \
    ../src/math/JonesMatrix.cpp \
    ../src/math/Mueller.cpp \
	../src/math/PhysMtr.cpp \
	../src/common/common.cpp \
	../src/tracing/Tracing.cpp \
    ../src/tracing/TracingConvex.cpp \
    ../src/tracing/TracingConcave.cpp \
	../src/CalcTimer.cpp \
	../src/Tracer.cpp \
	../src/geometry/clipper/clipper.cpp \
	../src/geometry/clipper/BeamClipper.cpp \
	../src/geometry/geometry_lib.cpp \
	../src/PhisBeam.cpp \
    ../src/particle/TiltedHexagonal.cpp \
    ../src/particle/HexagonalAggregate.cpp \
    ../src/particle/CertainAggregate.cpp \
    ../src/geometry/Polygon.cpp \
    ../src/geometry/Facet.cpp \
    ../src/particle/BulletRosette.cpp \
    ../src/particle/Bullet.cpp \
    ../src/common/Matrix4x4.cpp \
    ../src/common/MullerMatrix.cpp


HEADERS += \
	../src/bigint/BigInteger.hh \
	../src/bigint/BigIntegerAlgorithms.hh \
	../src/bigint/BigIntegerLibrary.hh \
	../src/bigint/BigIntegerUtils.hh \
	../src/bigint/BigUnsigned.hh \
	../src/bigint/BigUnsignedInABase.hh \
	../src/bigint/NumberlikeArray.hh \

SOURCES += \
	../src/bigint/BigInteger.cc \
	../src/bigint/BigIntegerAlgorithms.cc \
	../src/bigint/BigIntegerUtils.cc \
	../src/bigint/BigUnsigned.cc \
	../src/bigint/BigUnsignedInABase.cc \
