INCLUDEPATH += \
    $$SRC/math \
    $$SRC/common \
    $$SRC/particle \
    $$SRC/geometry \
<<<<<<< HEAD
#    $$SRC/scattering \
#    $$SRC/tracer \
#    $$SRC/incidence \
=======
    $$SRC/scattering \
    $$SRC/tracer \
    $$SRC/incidence \
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46
    $$SRC/bigint

SOURCES += \
    $$SRC/math/*.cpp \
<<<<<<< HEAD
    $$SRC/particle/*.cpp \
    $$SRC/geometry/*.cpp \
    $$SRC/common/*.cpp \
#    $$SRC/scattering/*.cpp \
#    $$SRC/incidence/*.cpp \
#    $$SRC/tracer/*.cpp \
=======
    $$SRC/handler/*.cpp \
    $$SRC/particle/*.cpp \
    $$SRC/geometry/*.cpp \
    $$SRC/common/*.cpp \
    $$SRC/scattering/*.cpp \
    $$SRC/incidence/*.cpp \
    $$SRC/tracer/*.cpp \
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46
    $$SRC/bigint/*.cc

HEADERS += \
    $$SRC/*.h \
<<<<<<< HEAD
=======
    $$SRC/handler/*.h \
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46
    $$SRC/math/*.hpp \
    $$SRC/math/*.h \
    $$SRC/particle/*.h \
    $$SRC/geometry/*.h \
    $$SRC/common/*.h \
<<<<<<< HEAD
#    $$SRC/scattering/*.h \
#    $$SRC/incidence/*.h \
#    $$SRC/tracer/*.h \
=======
    $$SRC/scattering/*.h \
    $$SRC/incidence/*.h \
    $$SRC/tracer/*.h \
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46
    $$SRC/bigint/*.hh
