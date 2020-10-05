TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
        eno_advection.cpp \
        grid2d.cpp \
        main.cpp

QMAKE_CXXFLAGS += -Xpreprocessor -fopenmp -lomp -I/usr/local/include

QMAKE_LFLAGS += -lomp

LIBS += -L /usr/local/lib /usr/local/lib/libomp.dylib


HEADERS += \
    cf_2.h \
    eno_advection.h \
    grid2d.h
