lib_LIBRARIES = libquantumdot.a
libquantumdot_a_SOURCES =  \
	QdotFci.cc \
	QdotFci.hpp \
	QdotInteraction.cc \
	QdotInteraction.hpp \
	QdotMatrixElements.hpp \
	QdotHilbertSpace.cc \
	QdotHilbertSpace.hpp \
	RadialPotential.hpp \
	RadialPotential.cc
LDADD = ../manybody/libmanybody.a ../ConfigFile/libConfigFile.a
LIBS = $(ARPACK_FLAGS) $(LAPACK_FLAGS)  $(FLENS_FLAGS)
INCLUDES = -I$(top_srcdir)/src/manybody -I$(HOME)/usr/include
libquantumdot_a_includedir = $(includedir)/quantumdot
libquantumdot_a_include_HEADERS = QdotFci.hpp \
	QdotInteraction.hpp \
	QdotMatrixElements.hpp \
	QdotHilbertSpace.hpp \
	RadialPotential.hpp



