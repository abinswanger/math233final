TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

QMAKE_CXXFLAGS += -fopenmp
QMAKE_LFLAGS += -fopenmp
LIBS += -fopenmp

SOURCES += \
        ../../lib/cpp/cf_2.cpp \
        ../../lib/cpp/eno_advection.cpp \
        ../../lib/cpp/grid2d.cpp \
        ../../lib/cpp/math_tools.cpp \
        ../../lib/cpp/sl_method.cpp \
        main.cpp \
        referencemap.cpp

HEADERS += \
    ../../lib/h/cf_2.h \
    ../../lib/h/eno_advection.h \
    ../../lib/h/grid2d.h \
    ../../lib/h/math_tools.h \
    ../../lib/h/sl_method.h \
    referencemap.h
