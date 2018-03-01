TEMPLATE = app

DESTDIR = ../bin

QMAKE_CXXFLAGS += -std=gnu++11
QMAKE_CXXFLAGS += -march=corei7 -msse4.2

TEMPLATE = app

CONFIG(release, debug|release): {
	DEFINES += _NDEBUG
	TARGET = btm
}

CONFIG(debug,	debug|release): {
	DEFINES += _DEBUG
	TARGET = btm_d
}

INCLUDEPATH += \
	../src \
	../src/math \
	../src/common \
	../src/particle \
	../src/geometry \
	../src/scattering \
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
