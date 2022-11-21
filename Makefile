rootPath = .
include ${rootPath}/include.mk

libSources = impl/*.c
libHeaders = inc/*.h
libTests = tests/*.c

stTafDependencies =  ${LIBDEPENDS}

all: all_libs all_progs
all_libs: dependencies ${LIBDIR}/libstTaf.a

all_progs: all_libs ${BINDIR}/stTafTests ${BINDIR}/maf_to_taf ${BINDIR}/taf_to_maf ${BINDIR}/taf_norm ${BINDIR}/taf_add_gap_bases

dependencies:
	mkdir -p ${LIBDIR} ${BINDIR}
	cd submodules/sonLib && PKG_CONFIG_PATH=${CWD}/lib/pkgconfig:${PKG_CONFIG_PATH} ${MAKE}
	mkdir -p ${BINDIR} ${LIBDIR} ${INCLDIR}
	rm -rf submodules/sonLib/bin/*.dSYM
	ln -f submodules/sonLib/bin/[a-zA-Z]* ${BINDIR}
	ln -f submodules/sonLib/lib/*.a ${LIBDIR}
	ln -f submodules/sonLib/lib/sonLib.a ${LIBDIR}/libsonLib.a

${LIBDIR}/libstTaf.a : ${libSources} ${libHeaders}  ${stTafDependencies}
	${CC} ${CFLAGS} ${LDFLAGS} -c ${libSources}
	${AR} rc libstTaf.a *.o
	mv libstTaf.a ${LIBDIR}/

${BINDIR}/stTafTests : ${libTests} ${LIBDIR}/libstTaf.a ${stTafDependencies}
	${CC} ${CFLAGS} ${LDFLAGS} -o ${BINDIR}/stTafTests ${libTests} ${LIBDIR}/libstTaf.a ${LDLIBS}

${BINDIR}/maf_to_taf : maf_to_taf.c ${LIBDIR}/libstTaf.a ${stTafDependencies}
	${CC} ${CFLAGS} -o ${BINDIR}/maf_to_taf maf_to_taf.c ${LIBDIR}/libstTaf.a ${LDLIBS}

${BINDIR}/taf_to_maf : taf_to_maf.c ${LIBDIR}/libstTaf.a ${stTafDependencies}
	${CC} ${CFLAGS} -o ${BINDIR}/taf_to_maf taf_to_maf.c ${LIBDIR}/libstTaf.a ${LDLIBS}

${BINDIR}/taf_norm : taf_norm.c ${LIBDIR}/libstTaf.a ${stTafDependencies}
	${CC} ${CFLAGS} -o ${BINDIR}/taf_norm taf_norm.c ${LIBDIR}/libstTaf.a ${LDLIBS}

${BINDIR}/taf_add_gap_bases : taf_add_gap_bases.cpp ${LIBDIR}/libstTaf.a ${stTafDependencies}
	${CXX} ${CPPFLAGS} ${CXXFLAGS} -o ${BINDIR}/taf_add_gap_bases taf_add_gap_bases.cpp ${LIBDIR}/libstTaf.a ${LDLIBS}

test : all
	${BINDIR}/stTafTests

clean :
	cd submodules/sonLib && ${MAKE} clean
	rm -rf *.o ${LIBDIR} ${BINDIR}

static :
	CFLAGS="$${CFLAGS} -static" \
	CPPFLAGS="$${CXXFLAGS} -static" \
	make all
