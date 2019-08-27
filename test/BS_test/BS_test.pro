QT += testlib
QT -= gui

CONFIG += qt console warn_on depend_includepath testcase c++11
CONFIG -= app_bundle

TEMPLATE = app

SOURCES +=  tst_po.cpp

QMAKE_CXXFLAGS += -std=gnu++11
QMAKE_CXXFLAGS += -march=corei7 -msse4.2

DEFINES += _TEST

CONFIG(release, debug|release): {
    TARGET = mbs
}

CONFIG(debug,	debug|release): {
    DEFINES += _DEBUG
    TARGET = mbs_d
}

INCLUDEPATH += \
    ../../src \
    ../../src/math \
    ../../src/handler \
    ../../src/common \
    ../../src/particle \
    ../../src/geometry \
    ../../src/scattering \
    ../../src/tracer \
    ../../src/bigint

SOURCES += \
    ../../src/Tracks.cpp \
    ../../src/math/*.cpp \
    ../../src/handler/*.cpp \
    ../../src/particle/*.cpp \
    ../../src/geometry/*.cpp \
    ../../src/common/*.cpp \
    ../../src/scattering/*.cpp \
    ../../src/bigint/*.cc

SOURCES -= ../../src/main.cpp
#message($$SOURCES)

HEADERS += \
    ../../src/*.h \
    ../../src/math/*.hpp \
    ../../src/math/*.h \
    ../../src/handler/*.h \
    ../../src/particle/*.h \
    ../../src/geometry/*.h \
    ../../src/common/*.h \
    ../../src/scattering/*.h \
    ../../src/bigint/*.hh

HEADERS -= ../../src/Beam.h

DISTFILES += \
    classes.qmodel \
    sequence.qmodel
