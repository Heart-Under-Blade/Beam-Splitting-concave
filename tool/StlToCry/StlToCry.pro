TEMPLATE = app
CONFIG += qt console arn_on depend_includepath c++11
CONFIG -= app_bundle
QT -= gui

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
        *.cpp
HEADERS += \
        *.h

SRC = ../../src

include(../../pro/MBS.pri)

SOURCES -= $$SRC/scattering/*.cpp $$SRC/incidence/*.cpp $$SRC/tracer/*.cpp ../VoronovParticleGenerator/*.cpp
HEADERS -= $$SRC/scattering/*.h $$SRC/incidence/*.h $$SRC/tracer/*.h

#message($$SOURCES)
