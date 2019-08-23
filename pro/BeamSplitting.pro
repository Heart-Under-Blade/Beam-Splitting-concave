TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle

DESTDIR = ../bin

VERSION = 2.0.0

QMAKE_CXXFLAGS += -std=gnu++11
QMAKE_CXXFLAGS += -march=corei7 -msse4.2

CONFIG(release, debug|release): {
    TARGET = mbs
}

CONFIG(debug,	debug|release): {
    DEFINES += _DEBUG
    TARGET = mbs_d
}

INCLUDEPATH += \
    ../src \
    ../src/math \
    ../src/incidence \
    ../src/common \
    ../src/particle \
    ../src/geometry \
    ../src/scattering \
    ../src/handler \
    ../src/tracer \
    ../src/bigint

SOURCES += \
    ../src/*.cpp \
    ../src/math/*.cpp \
    ../src/incidence/*.cpp \
    ../src/tracer/*.cpp \
    ../src/particle/*.cpp \
    ../src/geometry/*.cpp \
    ../src/handler/*.cpp \
    ../src/common/*.cpp \
    ../src/scattering/*.cpp \
    ../src/bigint/*.cc

SOURCES -= ../src/Beam.cpp

HEADERS += \
    ../src/*.h \
    ../src/math/*.hpp \
    ../src/incidence/*.h \
    ../src/tracer/*.h \
    ../src/math/*.h \
    ../src/particle/*.h \
    ../src/geometry/*.h \
    ../src/handler/*.h \
    ../src/common/*.h \
    ../src/scattering/*.h \
    ../src/bigint/*.hh

DISTFILES += \
    *.qmodel
