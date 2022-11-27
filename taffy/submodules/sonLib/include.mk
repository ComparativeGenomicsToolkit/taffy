# we do specific stuff for specific host for now.
HOSTNAME = $(shell hostname)
MACH = $(shell uname -m)
SYS =  $(shell uname -s)
PYTHON = python3

BINDIR ?= ${rootPath}/../sonLib/bin
LIBDIR ?= ${rootPath}/../sonLib/lib

##
# C compiler flags
#
# - enable PIC as some newer GCC versions are configured to link PIE by default,
#   which will not link with non-PIC objcts
CFLAGS += -fPIC -std=c99
CXXFLAGS += -fPIC

##
# CGL `standard' inc/impl
CPPFLAGS += -Iinc -Iimpl

##
# On ARM/clang (but not on Mac?) char is unsigned by default
# This breaks file parsing in sonLib which compares char to EOF (-1)
# https://stackoverflow.com/questions/19555248/clang-make-char-signed-by-default-on-arm
CFLAGS += -fsigned-char
CPPFLAGS += -fsigned-char

##
# For GCC, the C++ 11 aplication binary interface must match the version
# that HDF5 uses.  Change the ABI version can address errors about undefined functions that vary
# by string type, such as 
#   std::__cxx11::basic_string vs std::basic_string
#
# Specify one of:
#   CXX_ABI_DEF = -D_GLIBCXX_USE_CXX11_ABI=0
# or
#   CXX_ABI_DEF = -D_GLIBCXX_USE_CXX11_ABI=1
#
# in include.local.mk or the environment will change the API.
# You must do a make clean if you change this variable.

# special handling to get C++ ABI right on UCSC Centos 6 servers
ifeq (${CXX_ABI_DEF},)
ifneq ($(wildcard /etc/redhat-release),)
ifeq ($(shell hostname -d), gi.ucsc.edu)
    export CXX_ABI_DEF = -D_GLIBCXX_USE_CXX11_ABI=0
endif
endif
endif
CXXFLAGS += ${CXX_ABI_DEF}


# C++ compiler.
ifndef CXX
  ifeq ($(SYS),Darwin) #This is to deal with the Mavericks replacing gcc with clang fully
    CXX = clang++
  else
    CXX = g++
  endif
endif

# Compiler flags.

#Release compiler flags
CFLAGS_opt = -O3 -g -Wall --pedantic -funroll-loops -DNDEBUG 
CXXFLAGS_opt = -O3 -g -Wall -funroll-loops -DNDEBUG

#Debug flags (slow)
CFLAGS_dbg = -Wall -O0 -Werror --pedantic -g -fno-inline -UNDEBUG -Wno-error=unused-result
CXXFLAGS_dbg = -Wall -g -O0 -fno-inline -UNDEBUG

#Ultra Debug flags (really slow, checks for memory issues)
CFLAGS_ultraDbg = -Wall -Werror --pedantic -g -O1 -fno-inline -fno-omit-frame-pointer -fsanitize=address
CXXFLAGS_ultraDbg = -g -O1 -fno-inline -fno-omit-frame-pointer -fsanitize=address

#Profile flags
CFLAGS_prof = -Wall -Werror --pedantic -pg -O3 -g -Wno-error=unused-result
CXXFLAGS_prof = -pg -O3 -g -Wall -funroll-loops -DNDEBUG

#Flags to use
ifneq (${CGL_PROF},)
  CXXFLAGS += ${CXXFLAGS_prof}
  CFLAGS += ${CFLAGS_prof}
else ifneq (${CGL_DEBUG},)
  ifeq (${CGL_DEBUG},ultra)
    CXXFLAGS += ${CXXFLAGS_ultraDbg}
    CFLAGS += ${CFLAGS_ultraDbg}
  else
    CXXFLAGS += ${CXXFLAGS_dbg}
    CFLAGS += ${CFLAGS_dbg}
  endif
else
  CXXFLAGS += ${CXXFLAGS_opt}
  CFLAGS += ${CFLAGS_opt}
endif

# other commands
RANLIB ?= ranlib

# location of Tokyo cabinet
ifndef tokyoCabinetLib
HAVE_TOKYO_CABINET = $(shell pkg-config --exists tokyocabinet; echo $$?)
ifneq ($(wildcard /opt/local/include/tcbdb.h),)
   # OS/X with TC installed from MacPorts
   tcPrefix = /opt/local
   tokyoCabinetIncl = -I${tcPrefix}/include -DHAVE_TOKYO_CABINET=1
   tokyoCabinetLib = -L${tcPrefix}/lib -Wl,-rpath,${tcPrefix}/lib -ltokyocabinet -lz -lbz2 -lpthread -lm
else ifneq ($(wildcard /usr/local/include/tcbdb.h),)
   # /usr/local install (FreeBSD, etc)
   tcPrefix = /usr/local
   tokyoCabinetIncl = -I${tcPrefix}/include -DHAVE_TOKYO_CABINET=1
   tokyoCabinetLib = -L${tcPrefix}/lib -Wl,-rpath,${tcPrefix}/lib -ltokyocabinet -lz -lbz2 -lpthread -lm
else ifneq ($(wildcard /usr/include/tcbdb.h),)
   # /usr install (Ubuntu, and probably most Debain-based systems)
   tcPrefix = /usr
   tokyoCabinetIncl = -I${tcPrefix}/include -DHAVE_TOKYO_CABINET=1
   tokyoCabinetLib = -L${tcPrefix}/lib -Wl,-rpath,${tcPrefix}/lib -ltokyocabinet -lz -lbz2 -lpthread -lm
else ifeq (${HAVE_TOKYO_CABINET},0)
   # Install registered with pkg-config
   tokyoCabinetIncl = $(shell pkg-config --cflags tokyocabinet) -DHAVE_TOKYO_CABINET=1
   tokyoCabinetLib = $(shell pkg-config --libs-only-L tokyocabinet) -Wl,-rpath,$(shell pkg-config --variable=libdir tokyocabinet) $(shell pkg-config --libs-only-l --static tokyocabinet)
endif
endif

# location of Kyoto Tycoon
ifndef kyotoTycoonLib
HAVE_KYOTO_TYCOON = $(shell pkg-config --exists kyototycoon; echo $$?)
ifneq ($(wildcard /opt/local/include/ktcommon.h),)
   # OS/X with TC installed from MacPorts
   ttPrefix = /opt/local
   kyotoTycoonIncl = -I${ttPrefix}/include -DHAVE_KYOTO_TYCOON=1 
   kyotoTycoonLib = -L${ttPrefix}/lib -Wl,-rpath,${ttPrefix}/lib -lkyototycoon -lkyotocabinet -lz -lbz2 -lpthread -lm -lstdc++ 
else ifneq ($(wildcard /usr/local/include/ktcommon.h),)
   # /usr/local install (FreeBSD, etc)
   ttPrefix = /usr/local
   kyotoTycoonIncl = -I${ttPrefix}/include -DHAVE_KYOTO_TYCOON=1 
   kyotoTycoonLib = -L${ttPrefix}/lib -Wl,-rpath,${ttPrefix}/lib -lkyototycoon -lkyotocabinet -lz -lbz2 -lpthread -lm -lstdc++
else ifneq ($(wildcard /usr/include/ktcommon.h),)
   # /usr install (Ubuntu)
   ttPrefix = /usr
   kyotoTycoonIncl = -I${ttPrefix}/include -DHAVE_KYOTO_TYCOON=1 
   kyotoTycoonLib = -L${ttPrefix}/lib -Wl,-rpath,${ttPrefix}/lib -lkyototycoon -lkyotocabinet -lz -lbz2 -lpthread -lm -lstdc++
else ifeq (${HAVE_KYOTO_TYCOON},0)
   # Install registered with pkg-config
   kyotoTycoonIncl = $(shell pkg-config --cflags kyototycoon) -DHAVE_KYOTO_TYCOON=1
   kyotoTycoonLib = $(shell pkg-config --libs-only-L kyototycoon) -Wl,-rpath,$(shell pkg-config --variable=libdir kyototycoon) $(shell pkg-config --libs-only-l --static kyototycoon)
endif
endif

# location of hiredis
ifndef hiRedisLib
  HAVE_REDIS = $(shell pkg-config --exists hiredis; echo $$?)
  ifeq (${HAVE_REDIS},0)
    hiRedisLib = $(shell pkg-config --libs hiredis)
    incs = $(shell pkg-config --cflags hiredis)
    ifeq ($(findstring -I,${incs}),)
      # Broken 14.04 package
      hiRedisIncl = ${incs} -DHAVE_REDIS=1 -I/usr/include/hiredis
    else
      hiRedisIncl = ${incs} -DHAVE_REDIS=1
    endif
  endif
endif

dblibs = ${tokyoCabinetLib} ${kyotoTycoonLib} ${hiRedisLib} -lz -lm

