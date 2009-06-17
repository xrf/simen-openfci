bin_PROGRAMS = pairing slater_demo tabulate radial
pairing_SOURCES = pairing.cc
tabulate_SOURCES = tabulate.cc
radial_SOURCES = radial.cc
slater_demo_SOURCES = slater_demo.cc
LDADD = ../quantumdot/libquantumdot.a ../manybody/libmanybody.a ../ConfigFile/libConfigFile.a 
#LIBS = $(ARPACK_FLAGS) $(LAPACK_FLAGS)  $(FLENS_FLAGS)
INCLUDES += -I$(top_srcdir)/src/ConfigFile -I$(top_srcdir)/src/manybody -I$(top_srcdir)/src/quantumdot -I$(HOME)/usr/include

EXTRA_DIST = pairing.ini tabulate.ini
