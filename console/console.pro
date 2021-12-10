CONFIG += console
CONFIG += c++17
CONFIG -= app_bundle
CONFIG -= qt

QMAKE_CXXFLAGS_RELEASE += -O3
QMAKE_CXXFLAGS_RELEASE += -ffast-math

SOURCES += \
        ../Simulation/analysis.cpp \
        main.cpp

HEADERS += \
    ../Simulation/analysis.hpp \
    ../Simulation/node_2d.hpp \
    ../Simulation/node_3d.hpp \
    ../Simulation/node_base.hpp \
    ../Simulation/parameters.hpp \
    ../Simulation/random_thijs.hpp \
    ../Simulation/rndutils.hpp \
    ../Simulation/simulation.hpp \
    ../Simulation/voronoi.hpp \
    config_parser.h

DISTFILES += \
    config.ini
