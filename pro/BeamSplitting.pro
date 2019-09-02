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

CONFIG(debug, debug|release): {
    DEFINES += _DEBUG
    TARGET = mbs_d
}

DISTFILES += \
    *.qmodel

# Insert string below to include source files in your project
SRC = ../src

SOURCES += \
    ../src/main.cpp

include(MBS.pri)

#message($$SOURCES)
