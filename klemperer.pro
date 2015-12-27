#-------------------------------------------------
#
# Project created by QtCreator 2015-03-27T17:53:10
#
#-------------------------------------------------

CONFIG += c++11
QT       += core gui
QMAKE_CXXFLAGS_RELEASE += -march=native -Ofast

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = klemperer2
TEMPLATE = app


SOURCES += main.cpp\
        mainwindow.cpp \
    integrator.cpp \
    displaywidget.cpp

HEADERS  += mainwindow.h \
    integrator.h \
    displaywidget.h

FORMS    += mainwindow.ui

DISTFILES +=

RESOURCES += \
    shaders.qrc
