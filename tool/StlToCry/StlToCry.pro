TEMPLATE = app
CONFIG += qt console arn_on depend_includepath c++11
CONFIG -= app_bundle
QT -= gui

<<<<<<< HEAD
<<<<<<< HEAD
=======
DESTDIR = ../build

>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46
=======
DESTDIR = ../build

=======
>>>>>>> origin/feature/voronoi
>>>>>>> origin/refactor
QMAKE_CXXFLAGS += -std=gnu++11
QMAKE_CXXFLAGS += -march=corei7 -msse4.2

CONFIG(release, debug|release): {
    TARGET = stltocry
}

CONFIG(debug,	debug|release): {
    DEFINES += _DEBUG
    TARGET = stltocry_d
}

SOURCES += \
<<<<<<< HEAD
        *.cpp
HEADERS += \
        *.h
=======
        main.cpp Converter.cpp
HEADERS += \
        Converter.h
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46

SRC = ../../src

include(../../pro/MBS.pri)

<<<<<<< HEAD
SOURCES -= $$SRC/scattering/*.cpp $$SRC/incidence/*.cpp $$SRC/tracer/*.cpp ../VoronovParticleGenerator/*.cpp
HEADERS -= $$SRC/scattering/*.h $$SRC/incidence/*.h $$SRC/tracer/*.h

#message($$SOURCES)
=======
SOURCES -= $$SRC/scattering/*.cpp $$SRC/incidence/*.cpp $$SRC/tracer/*.cpp
HEADERS -= $$SRC/scattering/*.h $$SRC/incidence/*.h $$SRC/tracer/*.h

message($$SOURCES)
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46
