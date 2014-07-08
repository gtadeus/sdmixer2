#-------------------------------------------------
#
# Project created by QtCreator 2014-07-02T20:08:01
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = sdmixer2
TEMPLATE = app


SOURCES += main.cpp\
        sdmixer.cpp \
    pairfinder.cpp

HEADERS  += sdmixer.h \
    pairfinder.h

FORMS    += sdmixer.ui

LIBS += -lgsl -lgslcblas -lm
