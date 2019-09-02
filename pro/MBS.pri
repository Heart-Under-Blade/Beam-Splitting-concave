INCLUDEPATH += \
    $$SRC/math \
    $$SRC/common \
    $$SRC/handler \
    $$SRC/incidence \
    $$SRC/particle \
    $$SRC/geometry \
    $$SRC/scattering \
    $$SRC/tracer \
    $$SRC/bigint

SOURCES += \
    $$SRC/math/*.cpp \
    $$SRC/handler/*.cpp \
    $$SRC/incidence/*.cpp \
    $$SRC/particle/*.cpp \
    $$SRC/geometry/*.cpp \
    $$SRC/common/*.cpp \
    $$SRC/scattering/*.cpp \
    $$SRC/tracer/*.cpp \
    $$SRC/bigint/*.cc

HEADERS += \
    $$SRC/*.h \
    $$SRC/handler/*.h \
    $$SRC/incidence/*.h \
    $$SRC/math/*.hpp \
    $$SRC/math/*.h \
    $$SRC/particle/*.h \
    $$SRC/geometry/*.h \
    $$SRC/common/*.h \
    $$SRC/scattering/*.h \
    $$SRC/tracer/*.h \
    $$SRC/bigint/*.hh
