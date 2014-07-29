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


SOURCES += main.cpp\
        sdmixer.cpp \
    pairfinder.cpp \
    filter.cpp \
    reconstructor.cpp \
    settings.cpp

HEADERS  += sdmixer.h \
    pairfinder.h \
    filter.h \
    reconstructor.h \
    settings.h

FORMS    += sdmixer.ui

LIBS += -lm -llibtiff
#LIBS += -lgsl -lgslcblas -lm

#windows 8
win32:DEPENDPATH += C:\Progra~2\GnuWin32\
win32:DEPENDPATH += C:\Progra~2\GnuWin32\lib
win32:INCLUDEPATH += C:\Progra~2\GnuWin32\include

#windows 7
win32:DEPENDPATH += C:\Progra~1\GnuWin32\
win32:DEPENDPATH += C:\Progra~1\GnuWin32\lib
win32:INCLUDEPATH += C:\Progra~1\GnuWin32\include

win32:LIBS += -LC:\Progra~2\GnuWin32\lib
win32:LIBS += -LC:\Progra~2\GnuWin32\bin
win32:LIBS += -LC:\Progra~1\GnuWin32\lib
win32:LIBS += -LC:\Progra~1\GnuWin32\bin
win32:LIBS += -llibgsl
win32:LIBS += -llibgslcblas

 if(!debug_and_release|build_pass):CONFIG(debug, debug|release) {
    #mac:LIBS = $$member(LIBS, 0) $$member(LIBS, 1)_debug
    #win32:LIBS = $$member(LIBS, 0) $$member(LIBS, 1)d
 }
