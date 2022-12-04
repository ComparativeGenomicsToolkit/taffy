rootPath = .
include ${rootPath}/include.mk

srcDir = taffy/impl
libHeaders = taffy/inc/*.h
libTests = tests/*.c

all: all_libs all_progs
all_libs: ${LIBDIR}/libstTaf.a

all_progs: all_libs ${BINDIR}/taffy

sonLib: 
	mkdir -p ${LIBDIR} ${BINDIR}
	cd taffy/submodules/sonLib && PKG_CONFIG_PATH=${CWD}/lib/pkgconfig:${PKG_CONFIG_PATH} ${MAKE}
	mkdir -p ${BINDIR} ${LIBDIR} ${INCLDIR}
	rm -rf taffy/submodules/sonLib/bin/*.dSYM
	ln -f taffy/submodules/sonLib/lib/*.a ${LIBDIR}
	ln -f taffy/submodules/sonLib/lib/sonLib.a ${LIBDIR}/libsonLib.a

stTafDependencies = ${sonLibDir}/sonLib.a ${sonLibDir}/cuTest.a

${sonLibDir}/sonLib.a : sonLib

${sonLibDir}/cuTest.a : sonLib

${LIBDIR}/libstTaf.a : ${srcDir}/alignment_block.o ${srcDir}/line_iterator.o ${srcDir}/maf.o ${srcDir}/ond.o ${srcDir}/taf.o ${srcDir}/tai.o ${libHeaders} ${stTafDependencies}
	${AR} rc libstTaf.a ${srcDir}/alignment_block.o ${srcDir}/line_iterator.o ${srcDir}/maf.o ${srcDir}/ond.o ${srcDir}/taf.o ${srcDir}/tai.o
	mv libstTaf.a ${LIBDIR}/

${srcDir}/alignment_block.o : ${srcDir}/alignment_block.c ${libHeaders}
	${CC} ${CFLAGS} ${LDFLAGS} -o ${srcDir}/alignment_block.o -c ${srcDir}/alignment_block.c

${srcDir}/line_iterator.o : ${srcDir}/line_iterator.c ${libHeaders}
	${CC} ${CFLAGS} ${LDFLAGS} -o ${srcDir}/line_iterator.o -c ${srcDir}/line_iterator.c

${srcDir}/maf.o : ${srcDir}/maf.c ${libHeaders}
	${CC} ${CFLAGS} ${LDFLAGS} -o ${srcDir}/maf.o -c ${srcDir}/maf.c

${srcDir}/ond.o : ${srcDir}/ond.c ${libHeaders}
	${CC} ${CFLAGS} ${LDFLAGS} -o ${srcDir}/ond.o -c ${srcDir}/ond.c

${srcDir}/taf.o : ${srcDir}/taf.c ${libHeaders}
	${CC} ${CFLAGS} ${LDFLAGS} -o ${srcDir}/taf.o -c ${srcDir}/taf.c

${srcDir}/tai.o : ${srcDir}/tai.c ${libHeaders}
	${CC} ${CFLAGS} ${LDFLAGS} -o ${srcDir}/tai.o -c ${srcDir}/tai.c

${BINDIR}/stTafTests : ${libTests} ${LIBDIR}/libstTaf.a ${stTafDependencies}
	${CC} ${CFLAGS} ${LDFLAGS} -o ${BINDIR}/stTafTests ${libTests} ${LIBDIR}/libstTaf.a ${LDLIBS}

${BINDIR}/taffy : taf_norm.o taf_add_gap_bases.o taf_index.o taf_view.o taffy_main.o ${LIBDIR}/libstTaf.a ${libHeaders} ${stTafDependencies}
	${CXX} ${CPPFLAGS} ${CXXFLAGS} taf_norm.o taf_add_gap_bases.o taf_index.o taf_view.o taffy_main.o -o ${BINDIR}/taffy ${LIBDIR}/libstTaf.a ${LDLIBS}

taffy_main.o : taffy_main.cpp ${stTafDependencies} ${libHeaders}
	${CXX} ${CPPFLAGS} ${CXXFLAGS} -o taffy_main.o -c taffy_main.cpp

taf_norm.o : taf_norm.c ${stTafDependencies} ${libHeaders}
	${CC} ${CFLAGS} -o taf_norm.o -c taf_norm.c

taf_add_gap_bases.o : taf_add_gap_bases.cpp ${stTafDependencies} ${libHeaders}
	${CXX} ${CPPFLAGS} ${CXXFLAGS} -o taf_add_gap_bases.o -c taf_add_gap_bases.cpp

taf_index.o : taf_index.c ${stTafDependencies} ${libHeaders}
	${CC} ${CFLAGS} ${CFLAGS} -o taf_index.o -c taf_index.c

taf_view.o : taf_view.c ${stTafDependencies} ${libHeaders}
	${CC} ${CFLAGS} ${CFLAGS} -o taf_view.o -c taf_view.c

test : all ${BINDIR}/stTafTests
	${BINDIR}/stTafTests
	tests/tai/test_tai.py

clean :
	cd taffy/submodules/sonLib && ${MAKE} clean
	rm -rf *.o taffy/impl/*.o ${LIBDIR} ${BINDIR}

static :
	CFLAGS="$${CFLAGS} -static -march=nehalem" \
	CPPFLAGS="$${CXXFLAGS} -static -march=nehalem" \
	TAF_STATIC=1 \
	make all
