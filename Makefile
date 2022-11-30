rootPath = .
include ${rootPath}/include.mk

libSources = impl/*.c
libHeaders = inc/*.h
libTests = tests/*.c

stTafDependencies =  ${LIBDEPENDS}

all: all_libs all_progs
all_libs: dependencies ${LIBDIR}/stTaf.a

all_progs: all_libs ${BINDIR}/stTafTests ${BINDIR}/maf_to_taf ${BINDIR}/taf_to_maf ${BINDIR}/taf_norm ${BINDIR}/taf_add_gap_bases ${BINDIR}/taf_index ${BINDIR}/taf_find

dependencies:
	mkdir -p ${LIBDIR} ${BINDIR}
	cd submodules/sonLib && PKG_CONFIG_PATH=${CWD}/lib/pkgconfig:${PKG_CONFIG_PATH} ${MAKE}
	mkdir -p ${BINDIR} ${LIBDIR} ${INCLDIR}
	rm -rf submodules/sonLib/bin/*.dSYM
	ln -f submodules/sonLib/bin/[a-zA-Z]* ${BINDIR}
	ln -f submodules/sonLib/lib/*.a ${LIBDIR}

${LIBDIR}/stTaf.a : ${libSources} ${libHeaders}  ${stTafDependencies}
	${CC} ${CFLAGS} ${LDFLAGS} -c ${libSources}
	${AR} rc stTaf.a *.o
	mv stTaf.a ${LIBDIR}/

${BINDIR}/stTafTests : ${libTests} ${LIBDIR}/stTaf.a ${stTafDependencies}
	${CC} ${CFLAGS} ${LDFLAGS} -o ${BINDIR}/stTafTests ${libTests} ${LIBDIR}/stTaf.a ${LDLIBS}

${BINDIR}/maf_to_taf : maf_to_taf.c ${LIBDIR}/stTaf.a ${stTafDependencies}
	${CC} ${CFLAGS} -o ${BINDIR}/maf_to_taf maf_to_taf.c ${LIBDIR}/stTaf.a ${LDLIBS}

${BINDIR}/taf_to_maf : taf_to_maf.c ${LIBDIR}/stTaf.a ${stTafDependencies}
	${CC} ${CFLAGS} -o ${BINDIR}/taf_to_maf taf_to_maf.c ${LIBDIR}/stTaf.a ${LDLIBS}

${BINDIR}/taf_norm : taf_norm.c ${LIBDIR}/stTaf.a ${stTafDependencies}
	${CC} ${CFLAGS} -o ${BINDIR}/taf_norm taf_norm.c ${LIBDIR}/stTaf.a ${LDLIBS}

${BINDIR}/taf_add_gap_bases : taf_add_gap_bases.cpp ${LIBDIR}/stTaf.a ${stTafDependencies}
	${CXX} ${CPPFLAGS} ${CXXFLAGS} -o ${BINDIR}/taf_add_gap_bases taf_add_gap_bases.cpp ${LIBDIR}/stTaf.a ${LDLIBS}

${BINDIR}/taf_index : taf_index_main.c ${LIBDIR}/stTaf.a ${stTafDependencies}
	${CC} ${CFLAGS} ${CFLAGS} -o ${BINDIR}/taf_index taf_index_main.c ${LIBDIR}/stTaf.a ${LDLIBS}

${BINDIR}/taf_find : taf_find_main.c ${LIBDIR}/stTaf.a ${stTafDependencies}
	${CC} ${CFLAGS} ${CFLAGS} -o ${BINDIR}/taf_find taf_find_main.c ${LIBDIR}/stTaf.a ${LDLIBS}

test : all
	${BINDIR}/stTafTests

clean :
	cd submodules/sonLib && ${MAKE} clean
	rm -rf *.o ${LIBDIR} ${BINDIR}

static :
	CFLAGS="$${CFLAGS} -static" \
	CPPFLAGS="$${CXXFLAGS} -static" \
	make all
