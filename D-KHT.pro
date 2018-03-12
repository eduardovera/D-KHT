QT -= core
QT -= gui
QT += widgets

TARGET = D-KHT
CONFIG += console
CONFIG += c++11
CONFIG -= app_bundle
CONFIG += qt

TEMPLATE = app

DEFINES += DLIB_PNG_SUPPORT=1

INCLUDEPATH += include

win32:CONFIG(release, debug|release): LIBS += -L$$PWD/libs/dlib/lib/ -ldlib
else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/libs/dlib/lib/ -ldlibd

unix:!macx: LIBS += -ldlib

INCLUDEPATH += $$PWD/libs/dlib
DEPENDPATH += $$PWD/libs/dlib


SOURCES += main.cpp \
    accumulatorball_t.cpp \
    hough.cpp \
    voting.cpp \
    quadtree_t.cpp \
    logger.cpp

HEADERS += \
    accum_ball_cell_t.h \
    accum_cell_t.h \
    accumulatorball_t.h \
    bin_t.h \
    hough.h \
    kernel_t.h \
    peak_detection.h \
    plane_t.h \
    reader_file.h \
    settings.h \
    voting.h \
    quadtree_t.h \
    sat.h \
    logger.h

DISTFILES +=

