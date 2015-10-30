#
AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS} 

lib_LTLIBRARIES = libopenrelay-@OPENRELAY_API_VERSION@.la
libopenrelay_@OPENRELAY_API_VERSION@_la_SOURCES = filter.cpp relay.cpp types.cpp \
acquisition/acquisition.cpp acquisition/encoder.cpp \
control/control.cpp \
display/nscreen.cpp \
measures/measures.cpp  measures/rna.cpp measures/WangSunFrequency.cpp measures/fuzzy.cpp measures/wavelet.cpp measures/switch.cpp\
logic/logic.cpp \
protection/protection.cpp \
protection/generator.cpp \
savedata/oscillography.cpp\
acquisition/IEC61850.cpp
#interface/openrelayxml.cpp 
libopenrelay_@OPENRELAY_API_VERSION@_la_LDFLAGS = -version-info $(OPENRELAY_SO_VERSION)
openrelay_includedir = $(includedir)/openrelay-$(OPENRELAY_API_VERSION)
nobase_openrelay_include_HEADERS = channel.h filter.h parameters.h relay.h types.h \
acquisition/acquisition.h acquisition/encoder.h \
control/control.h \
display/nscreen.h \
logic/logic.h \
measures/measures.h measures/rna.h measures/WangSunFrequency.h measures/fuzzy.h measures/wavelet.h measures/switch.h \
protection/protection.h \
protection/generator.h \
savedata/oscillography.h\
acquisition/IEC61850.h
#interface/openrelayxml.h

openrelay_libincludedir = $(libdir)/openrelay-$(OPENRELAY_API_VERSION)/include

#bin_PROGRAMS = openrelay

#openrelay_SOURCES = interface/openrelay.cppcd i	

# openrelay_CPPFLAGS = -I../orelay  ${CPPFLAGS}
#openrelay_CPPFLAGS = -shared  ${CPPFLAGS}
#openrelay_LDFLAGS = -shared  ${LDFLAGS}    
#openrelay_LDADD = ./libopenrelay-${OPENRELAY_API_VERSION}.la -lm -lncurses -lfann -lpthread -lboost_thread 









