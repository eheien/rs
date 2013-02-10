TARGET   = siab
QT += opengl
MACPORTBASE    = /opt/local
INCLUDEPATH = $${MACPORTBASE}/include/qwt
INCLUDEPATH += $${MACPORTBASE}/include
LIBS        += -L$${MACPORTBASE}/lib -lqwt -lsundials_cvode -lsundials_nvecserial

HEADERS = \
    RateStateSimWindow.h \
    RateState.h \

SOURCES = \
    RateState.cpp \
    RateStateSimWindow.cpp \
    main.cpp

