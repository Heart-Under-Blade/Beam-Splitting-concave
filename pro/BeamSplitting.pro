TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
#CONFIG -= qt

DESTDIR = ../bin

VERSION = 1.2.0

QMAKE_CXXFLAGS += -std=gnu++11
QMAKE_CXXFLAGS += -march=corei7 -msse4.2

CONFIG(release, debug|release): {
	DEFINES += _NDEBUG
    TARGET = mbs
}

CONFIG(debug,	debug|release): {
	DEFINES += _DEBUG
    TARGET = mbs_d
}

INCLUDEPATH += \
	../src \
	../src/math \
	../src/common \
	../src/particle \
	../src/geometry \
	../src/scattering \
    ../src/tracer \
	../src/bigint

SOURCES += \
	../src/*.cpp \
	../src/math/*.cpp \
	../src/particle/*.cpp \
	../src/geometry/*.cpp \
	../src/common/*.cpp \
	../src/scattering/*.cpp \
    ../src/bigint/*.cc

HEADERS += \
	../src/*.h \
	../src/math/*.hpp \
	../src/math/*.h \
	../src/particle/*.h \
	../src/geometry/*.h \
	../src/common/*.h \
	../src/scattering/*.h \
    ../src/bigint/*.hh

DISTFILES += \
	classes.qmodel \
    sequence.qmodel
