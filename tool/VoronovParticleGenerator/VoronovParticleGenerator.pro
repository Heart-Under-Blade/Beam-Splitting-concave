TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
QT -= gui

#DESTDIR = ../bin

QMAKE_CXXFLAGS += -std=gnu++11
QMAKE_CXXFLAGS += -march=corei7 -msse4.2

CONFIG(release, debug|release): {
    TARGET = vpg
}

CONFIG(debug, debug|release): {
    DEFINES += _DEBUG
    TARGET = vpg_d
}

INCLUDEPATH += ../../src

SRC = ../../src

SOURCES += \
    *.cpp

include(../../pro/MBS.pri)

HEADERS += *.h
