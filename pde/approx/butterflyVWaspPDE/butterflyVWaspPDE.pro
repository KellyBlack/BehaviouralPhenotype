TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt
LIBS += -lpthread
QMAKE_CXXFLAGS += -O3

SOURCES += \
    approximationbase.cpp \
        main.cpp \
    pdeSolver.cpp \
    butterflies.cpp \
    numericaltrials.cpp \
    rungakutta45.cpp

HEADERS += \
    approximationbase.h \
    legendre.h \
    limInf.h \
    mainRoutines.h \
    util.h \
    lu_decomp.h \
    pdeSolver.h \
    butterflies.h \
    numericaltrials.h \
    rungakutta45.h
