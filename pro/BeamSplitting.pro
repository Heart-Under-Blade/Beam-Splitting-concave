TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

DESTDIR = ../bin

QMAKE_CXXFLAGS += -std=gnu++11
QMAKE_CXXFLAGS += -march=corei7 -msse4.2

CONFIG(release, debug|release): DEFINES += _NDEBUG
CONFIG(debug, debug|release): DEFINES += _DEBUG

INCLUDEPATH += \
	../src \
	../src/math \
	../src/particle \
	../src/common \
	../src/tracing

HEADERS += \
	../src/Beam.h \
	../src/Intersection.h \
	../src/test.h \
    ../src/intrinsic/intrinsics.h \
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
    ../src/common/types.h \
    ../src/common/vector_lib.h \
	../src/common/macro.h \
    ../src/tracing/TracingConvex.h \
    ../src/tracing/TracingConcave.h \
	../src/CalcTimer.h \
	../src/Tracer.h \
    ../src/clipper.hpp

SOURCES += \
    ../src/Beam.cpp \
	../src/Intersection.cpp \
	../src/main.cpp \
    ../src/intrinsic/intrinsics.cpp \
    ../src/math/compl.cpp \
	../src/math/matrix.cpp \
    ../src/particle/Hexagonal.cpp \
	../src/particle/ConcaveHexagonal.cpp \
    ../src/particle/Particle.cpp \
    ../src/math/JonesMatrix.cpp \
    ../src/math/Mueller.cpp \
	../src/math/PhysMtr.cpp \
    ../src/common/vector_lib.cpp \
	../src/common/common.cpp \
	../src/tracing/Tracing.cpp \
    ../src/tracing/TracingConvex.cpp \
    ../src/tracing/TracingConcave.cpp \
	../src/CalcTimer.cpp \
	../src/Tracer.cpp \
    ../src/clipper.cpp

