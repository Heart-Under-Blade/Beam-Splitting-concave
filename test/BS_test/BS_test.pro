QT += testlib
QT -= gui

CONFIG += qt console warn_on depend_includepath testcase c++11
CONFIG -= app_bundle

TEMPLATE = app

QMAKE_CXXFLAGS += -std=gnu++11
QMAKE_CXXFLAGS += -march=corei7 -msse4.2

DEFINES += _TEST

CONFIG(release, debug|release): {
    TARGET = mbs_tst
}

CONFIG(debug,	debug|release): {
    DEFINES += _DEBUG
    TARGET = mbs_tst_d
}

SRC = ../../src

INCLUDEPATH += $$SRC

SOURCES += tst_po.cpp  \
    $$SRC/Tracks.cpp \

HEADERS += \
    $$SRC/*.h \

include(../../pro/MBS.pri)
