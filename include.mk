SHELL = /bin/bash

##
# Users can set CPPFLAGS, CFLAGS, LIBS to reference external packages.  These
# can be set in environment of config.local.mk.  LDLIBS should not be modified,
# as it is not seen my kyoto configure.
##

# special handling to get C++ ABI right on UCSC Centos 6 servers
ifeq (${CXX_ABI_DEF},)
ifneq ($(wildcard /etc/redhat-release),)
ifeq ($(shell hostname -d), gi.ucsc.edu)
    export CXX_ABI_DEF = -D_GLIBCXX_USE_CXX11_ABI=0
endif
endif
endif


#Location of sonLib
BINDIR = ${rootPath}/bin
LIBDIR = ${rootPath}/lib
INCLDIR = ${rootPath}/include

#Modify this variable to set the location of sonLib
sonLibRootDir ?= ${rootPath}/taffy/submodules/sonLib
sonLibDir=${sonLibRootDir}/lib

include ${sonLibRootDir}/include.mk

#Turn asserts back on in spite of sonLib
#https://github.com/ComparativeGenomicsToolkit/cactus/issues/235
CFLAGS += -UNDEBUG

ifndef TARGETOS
  TARGETOS := $(shell uname -s)
endif

# Hack to include openmp on os x after "brew install lomp
ifeq ($(TARGETOS), Darwin)
	CFLAGS+= -Xpreprocessor -fopenmp -lomp
	CXXFLAGS+= -Xpreprocessor -fopenmp -lomp
else
	CFLAGS+= -fopenmp
	CXXFLAGS+= -fopenmp
endif

# Hack in ARM support
# Toggle on if "arm" is set, or if uname -m returns aarch64
ifeq ($(shell uname -m || true), aarch64)
	arm=1
endif
ifeq ($(shell arch || true), aarch64)
	arm=1
endif
ifdef arm
# flags to build abpoa
export armv8 = 1
export aarch64 = 1
# flags to include simde abpoa in cactus on ARM
CFLAGS+= -march=armv8-a+simd
else
# flags to build abpoa
export avx2 = 1
endif
# flags needed to include simde abpoa in cactus on any architecture
CFLAGS+= -D__AVX2__ -DUSE_SIMDE -DSIMDE_ENABLE_NATIVE_ALIASES

dataSetsPath=/Users/benedictpaten/Dropbox/Documents/work/myPapers/genomeCactusPaper/dataSets

inclDirs = taffy/inc taffy/submodules/sonLib/C/inc taffy/submodules/sonLib/externalTools/cutest

# htslib will be used only if it can be found via pkg-config
# without it, taf won't support bgzipped input (but is otherwise the same)
HAVE_HTSLIB = $(shell pkg-config --exists htslib; echo $$?)
ifeq (${HAVE_HTSLIB},0)
	HTSLIB_CFLAGS = $(shell pkg-config htslib --cflags) -DUSE_HTSLIB
	ifdef TAF_STATIC
		HTSLIB_LIBS = $(shell pkg-config htslib --libs --static) -llzma
	else
		HTSLIB_LIBS = $(shell pkg-config htslib --libs)
	endif
endif

CFLAGS += ${inclDirs:%=-I${rootPath}/%} -I${LIBDIR} -I${rootPath}/include  ${HTSLIB_CFLAGS}
CXXFLAGS += ${inclDirs:%=-I${rootPath}/%} -I${LIBDIR} -I${rootPath}/include ${HTSLIB_CFLAGS}

LDLIBS += ${HTSLIB_LIBS}

# libraries can't be added until they are build, so add as to LDLIBS until needed
sonLibLibs = ${sonLibDir}/sonLib.a ${sonLibDir}/cuTest.a

# optional hal support toggled on by setting HALDIR (ex to ../hal)
ifdef HALDIR
	LDLIBS += ${HALDIR}/lib/libHalBlockViz.a ${HALDIR}/lib/libHalLiftover.a ${HALDIR}/lib/libHalLod.a ${HALDIR}/lib/libHalMaf.a ${HALDIR}/lib/libHal.a
	CFLAGS += -I${HALDIR}/api/inc -I${HALDIR}/blockViz/inc -DUSE_HAL
	CXXFLAGS += -I${HALDIR}/api/inc -I${HALDIR}/blockViz/inc -DUSE_HAL
	CXX = h5c++

# This bit copied from HAL: todo -- we should probably just use hal's include.mk directly
# add compiler flag and kent paths if udc is enabled
# relies on KENTSRC containing path to top level kent/src dir
# and MACHTYPE being specified.
# This MUST follow PHAST defs, as they both have a gff.h
	ifdef ENABLE_UDC
#  Find htslib as in kent/src/inc/common.mk:
		MACHTYPE = x86_64
		CXXFLAGS += -DENABLE_UDC-I${KENTSRC}/inc -I${KENTSRC}/htslib -pthread
		LDLIBS += ${KENTSRC}/lib/${MACHTYPE}/jkweb.a ${KENTSRC}/htslib/libhts.a -lcurl -lssl -lcrypto -pthread
	endif

endif

# note: the CACTUS_STATIC_LINK_FLAGS below can generally be empty -- it's used by the static builder script only
LDLIBS += ${sonLibLibs} ${LIBS} -L${rootPath}/lib -Wl,-rpath,${rootPath}/lib -lz -lbz2 -lpthread -lm -lstdc++ -lm ${CACTUS_STATIC_LINK_FLAGS}

