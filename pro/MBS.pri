INCLUDEPATH += \
    $$SRC/math \
    $$SRC/common \
    $$SRC/particle \
    $$SRC/geometry \
    $$SRC/scattering \
    $$SRC/tracer \
    $$SRC/incidence \
    $$SRC/bigint

SOURCES += \
    $$SRC/math/*.cpp \
    $$SRC/particle/*.cpp \
    $$SRC/geometry/*.cpp \
    $$SRC/common/*.cpp \
    $$SRC/scattering/*.cpp \
    $$SRC/incidence/*.cpp \
    $$SRC/tracer/*.cpp \
    $$SRC/bigint/*.cc

HEADERS += \
    $$SRC/*.h \
    $$SRC/math/*.hpp \
    $$SRC/math/*.h \
    $$SRC/particle/*.h \
    $$SRC/geometry/*.h \
    $$SRC/common/*.h \
    $$SRC/scattering/*.h \
    $$SRC/incidence/*.h \
    $$SRC/tracer/*.h \
    $$SRC/bigint/*.hh
