TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt
LIBS += -lpthread
QMAKE_CXXFLAGS += -O3

SOURCES += \
        main.cpp \
    pdeSolver.cpp \
    butterflies.cpp \
    numericaltrials.cpp

HEADERS += \
    legendre.h \
    util.h \
    lu_decomp.h \
    pdeSolver.h \
    butterflies.h \
    numericaltrials.h
