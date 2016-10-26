TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

DESTDIR = $$PWD

QMAKE_CXXFLAGS += -std=gnu++11
QMAKE_CXXFLAGS += -march=corei7 -msse4.2

INCLUDEPATH += \
	../src \
	../src/math \
	../src/particle \
	../src/common \
	../src/tracing

SOURCES += \
    main.cpp \
	../src/common/vector_lib.cpp \
    ../src/particle/ConcaveHexagonal.cpp \
    ../src/particle/Hexagonal.cpp \
    ../src/particle/Particle.cpp

HEADERS += \
	../src/common/types.h \
	../src/common/vector_lib.h \
    ../src/particle/ConcaveHexagonal.h \
    ../src/particle/Hexagonal.h \
    ../src/particle/Particle.h
