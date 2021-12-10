QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets printsupport

CONFIG += c++17

# You can also make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

SOURCES += \
    main.cpp \
    mainwindow.cpp \
    qcustomplot.cpp

HEADERS += \
    ../Simulation/node_2d.hpp \
    ../Simulation/node_3d.hpp \
    ../Simulation/node_base.hpp \
    ../Simulation/parameters.hpp \
    ../Simulation/random_thijs.hpp \
    ../Simulation/rndutils.hpp \
    ../Simulation/setup.hpp \
    ../Simulation/simulation.hpp \
    ../Simulation/voronoi.hpp \
    ../Simulation/voronoi_tools.hpp \
    Simulation/node_2d.hpp \
    Simulation/node_3d.hpp \
    Simulation/node_base.hpp \
    Simulation/parameters.hpp \
    Simulation/random_thijs.hpp \
    Simulation/rndutils.hpp \
    Simulation/setup.hpp \
    Simulation/simulation.hpp \
    Simulation/voronoi.hpp \
    Simulation/voronoi_tools.hpp \
    mainwindow.h \
    qcustomplot.h

FORMS += \
    mainwindow.ui

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target

CONFIG+=force_debug_info CONFIG+=separate_debug_info

CONFIG(release, debug|release) {
    CONFIG += optimize_full
}

QMAKE_CXXFLAGS_RELEASE += -ffast-math

TARGET = ovr
