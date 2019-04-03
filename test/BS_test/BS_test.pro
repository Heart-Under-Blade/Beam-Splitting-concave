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

SRC = ../../src

SOURCES +=
    $$SRC/*.cpp \

include(../../pro/MBS.pri)
