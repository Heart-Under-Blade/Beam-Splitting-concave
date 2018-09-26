QT += testlib
QT -= gui

CONFIG += qt console warn_on depend_includepath testcase c++11
CONFIG -= app_bundle

TEMPLATE = app

SOURCES += tst_po.cpp

QMAKE_CXXFLAGS += -std=gnu++11
QMAKE_CXXFLAGS += -march=corei7 -msse4.2

CONFIG(release, debug|release): {
    DEFINES += _NDEBUG
    TARGET = mbs_tst
}

CONFIG(debug,	debug|release): {
    DEFINES += _DEBUG
    TARGET = mbs_tst_d
}

INCLUDEPATH += \
    ../../src/math \
    ../../src/common \
    ../../src/particle \
    ../../src/geometry \
    ../../src/scattering \
    ../../src/tracer \
    ../../src/bigint \
    ../../src/incidence

SOURCES += \
    ../../src/math/*.cpp \
    ../../src/particle/*.cpp \
    ../../src/geometry/*.cpp \
    ../../src/common/*.cpp \
    ../../src/scattering/*.cpp \
    ../../src/incidence/*.cpp \
    ../../src/bigint/*.cc

message($$SOURCES)

HEADERS += \
    ../../src/*.h \
    ../../src/math/*.hpp \
    ../../src/math/*.h \
    ../../src/particle/*.h \
    ../../src/geometry/*.h \
    ../../src/common/*.h \
    ../../src/scattering/*.h \
    ../../src/bigint/*.hh

DISTFILES += \
    classes.qmodel \
    sequence.qmodel
