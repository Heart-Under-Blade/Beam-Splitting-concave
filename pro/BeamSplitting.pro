TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

DESTDIR = ../bin

QMAKE_CXXFLAGS += -std=gnu++11
QMAKE_CXXFLAGS += -march=corei7 -msse4.2


INCLUDEPATH += ../src/math

HEADERS += \
    ../src/Beam.h \
    ../src/global.h \
    ../src/Hexagon.h \
    ../src/Intersection.h \
    ../src/JonesMatrix.h \
    ../src/Mueller.hpp \
    ../src/Particle.h \
    ../src/PhysMtr.hpp \
    ../src/test.h \
    ../src/Tracing.h \
    ../src/types.h \
    ../src/vector_lib.h \
    ../src/intrinsic/intrinsics.h \
    ../src/math/compl.hpp \
    ../src/math/matrix.hpp \
    ../src/math/service.hpp

SOURCES += \
    ../src/Beam.cpp \
    ../src/Hexagon.cpp \
    ../src/Intersection.cpp \
    ../src/JonesMatrix.cpp \
    ../src/main.cpp \
    ../src/Mueller.cpp \
    ../src/Particle.cpp \
    ../src/PhysMtr.cpp \
    ../src/Tracing.cpp \
    ../src/vector_lib.cpp \
    ../src/intrinsic/intrinsics.cpp \
    ../src/math/compl.cpp \
    ../src/math/matrix.cpp

