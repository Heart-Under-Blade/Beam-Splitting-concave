TEMPLATE = app
CONFIG += qt console arn_on depend_includepath c++11
CONFIG -= app_bundle
QT -= gui

QMAKE_CXXFLAGS += -std=gnu++11
QMAKE_CXXFLAGS += -march=corei7 -msse4.2

CONFIG(release, debug|release): {
    DEFINES += _NDEBUG
    TARGET = stltocry
}

CONFIG(debug,	debug|release): {
    DEFINES += _DEBUG
    TARGET = stltocry_d
}

SOURCES += \
        main.cpp

SRC = ../../src

include(../../pro/MBS.pri)
