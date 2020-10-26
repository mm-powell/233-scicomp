TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
        grid2d.cpp \
        main.cpp \
        math_tools.cpp \
        sl_method.cpp

QMAKE_CXXFLAGS += -Xpreprocessor -fopenmp -lomp -I/usr/local/include

QMAKE_LFLAGS += -lomp

LIBS += -L /usr/local/lib /usr/local/lib/libomp.dylib

HEADERS += \
    cf_2.h \
    grid2d.h \
    math_tools.h \
    sl_method.h
