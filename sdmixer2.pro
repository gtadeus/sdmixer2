#-------------------------------------------------
#
# Project created by QtCreator 2014-07-02T20:08:01
#
#-------------------------------------------------

QT       += core gui
QT += xml

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = sdmixer2
TEMPLATE = app
CONFIG += c++11
win32:RC_FILE = run.rc

QMAKE_CFLAGS  -= -O2
QMAKE_CFLAGS  -= -O1
QMAKE_CXXFLAGS  -= -O2
QMAKE_CXXFLAGS  -= -O1
QMAKE_CFLAGS  =  -O3
QMAKE_LFLAGS  =  -O3
QMAKE_CXXFLAGS  = -O3


SOURCES += main.cpp\
sdmixer.cpp \
    pairfinder.cpp \
    filter.cpp \
    reconstructor.cpp \
    settings.cpp \
    mapped_file.cpp



HEADERS  += sdmixer.h \
    pairfinder.h \
    filter.h \
    reconstructor.h \
    settings.h \


FORMS    += sdmixer.ui

LIBS += -llibtiff


win32:INCLUDEPATH += C:\Build\include

win32:LIBS += -LC:\Build\libs



 if(!debug_and_release|build_pass):CONFIG(debug, debug|release) {
    #mac:LIBS = $$member(LIBS, 0) $$member(LIBS, 1)_debug
    #win32:LIBS = $$member(LIBS, 0) $$member(LIBS, 1)d
 }
