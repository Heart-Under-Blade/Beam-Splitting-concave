TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

DESTDIR = $$PWD

QMAKE_CXXFLAGS += -std=gnu++11
QMAKE_CXXFLAGS += -march=corei7 -msse4.2

INCLUDEPATH += \
	../../src \
	../../src/math \
	../../src/particle \
	../../src/common \
	../../src/tracing \
	../../src/geometry

SOURCES += \
	main.cpp \
	../../src/particle/ConcaveHexagonal.cpp \
	../../src/particle/Hexagonal.cpp \
	../../src/particle/HexagonalAggregate.cpp \
	../../src/particle/Particle.cpp \
	../../src/geometry/geometry_lib.cpp

HEADERS += \
	../../src/particle/ConcaveHexagonal.h \
	../../src/particle/Hexagonal.h \
	../../src/particle/HexagonalAggregate.cpp \
	../../src/particle/Particle.h \
	../../src/geometry/geometry_lib.h
