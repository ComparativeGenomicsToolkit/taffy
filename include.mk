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
sonLibRootDir ?= ${rootPath}/submodules/sonLib
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
else
	CFLAGS+= -fopenmp
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

inclDirs = inc submodules/sonLib/C/inc submodules/sonLib/externalTools/cutest

CPPFLAGS += ${inclDirs:%=-I${rootPath}/%} -I${LIBDIR} -I${rootPath}/include

# libraries can't be added until they are build, so add as to LDLIBS until needed
sonLibLibs = ${sonLibDir}/sonLib.a ${sonLibDir}/cuTest.a

# optional hal support toggled on by setting HALDIR (ex to ../hal)
ifdef HALDIR
	LDLIBS += ${HALDIR}/lib/libHalBlockViz.a ${HALDIR}/lib/libHalLiftover.a ${HALDIR}/lib/libHalLod.a ${HALDIR}/lib/libHalMaf.a ${HALDIR}/lib/libHal.a
	CFLAGS += -I${HALDIR}/api/inc -I${HALDIR}/blockViz/inc -DUSE_HAL
	CPPFLAGS += -I${HALDIR}/api/inc -I${HALDIR}/blockViz/inc -DUSE_HAL
	CXX = h5c++
endif

# note: the CACTUS_STATIC_LINK_FLAGS below can generally be empty -- it's used by the static builder script only
LDLIBS += ${sonLibLibs} ${LIBS} -L${rootPath}/lib -Wl,-rpath,${rootPath}/lib -lz -lbz2 -lpthread -lm -lstdc++ -lm ${CACTUS_STATIC_LINK_FLAGS}
LIBDEPENDS = ${sonLibDir}/sonLib.a ${sonLibDir}/cuTest.a

