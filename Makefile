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
	cd ${sonLibRootDir} && PKG_CONFIG_PATH=${CWD}/lib/pkgconfig:${PKG_CONFIG_PATH} ${MAKE}
	mkdir -p ${BINDIR} ${LIBDIR} ${INCLDIR}
	rm -rf ${sonLibRootDir}/bin/*.dSYM
	ln -f ${sonLibDir}/*.a ${LIBDIR}
	ln -f ${sonLibDir}/sonLib.a ${LIBDIR}/libsonLib.a

abPOA:
	mkdir -p ${LIBDIR} ${INCLDIR}
	cd taffy/submodules/abPOA && ${MAKE}
	ln -f taffy/submodules/abPOA/lib/*.a ${LIBDIR}
	ln -f taffy/submodules/abPOA/include/*.h ${INCLDIR}
	rm -fr ${INCLDIR}/simde && cp -r taffy/submodules/abPOA/include/simde ${INCLDIR}

${LIBDIR}/libabpoa.a : abPOA

${sonLibDir}/sonLib.a : sonLib

${sonLibDir}/cuTest.a : sonLib

stTafDependencies = ${sonLibDir}/sonLib.a ${sonLibDir}/cuTest.a ${LIBDIR}/libabpoa.a

${LIBDIR}/libstTaf.a : ${libTests} ${libHeaders} ${srcDir}/alignment_block.o ${srcDir}/line_iterator.o ${srcDir}/maf.o ${srcDir}/paf.o ${srcDir}/ond.o ${srcDir}/taf.o ${srcDir}/add_gap_bases.o ${srcDir}/merge_adjacent_alignments.o ${srcDir}/prefix_sort.o ${srcDir}/wiggle.o ${srcDir}/tai.o ${libHeaders} ${stTafDependencies}
	${AR} rc libstTaf.a ${srcDir}/alignment_block.o ${srcDir}/line_iterator.o ${srcDir}/maf.o ${srcDir}/paf.o ${srcDir}/ond.o ${srcDir}/taf.o ${srcDir}/add_gap_bases.o ${srcDir}/merge_adjacent_alignments.o ${srcDir}/prefix_sort.o ${srcDir}/wiggle.o ${srcDir}/tai.o
	mv libstTaf.a ${LIBDIR}/

${srcDir}/alignment_block.o : ${srcDir}/alignment_block.c ${libHeaders}
	${CC} ${CFLAGS} ${LDFLAGS} -o ${srcDir}/alignment_block.o -c ${srcDir}/alignment_block.c

${srcDir}/line_iterator.o : ${srcDir}/line_iterator.c ${libHeaders}
	${CC} ${CFLAGS} ${LDFLAGS} -o ${srcDir}/line_iterator.o -c ${srcDir}/line_iterator.c

${srcDir}/maf.o : ${srcDir}/maf.c ${libHeaders}
	${CC} ${CFLAGS} ${LDFLAGS} -o ${srcDir}/maf.o -c ${srcDir}/maf.c

${srcDir}/paf.o : ${srcDir}/paf.c ${libHeaders}
	${CC} ${CFLAGS} ${LDFLAGS} -o ${srcDir}/paf.o -c ${srcDir}/paf.c

${srcDir}/ond.o : ${srcDir}/ond.c ${libHeaders}
	${CC} ${CFLAGS} ${LDFLAGS} -o ${srcDir}/ond.o -c ${srcDir}/ond.c

${srcDir}/add_gap_bases.o : ${srcDir}/add_gap_bases.cpp ${libHeaders}
	${CXX} ${CPPFLAGS} ${CXXFLAGS} -o ${srcDir}/add_gap_bases.o -c ${srcDir}/add_gap_bases.cpp

${srcDir}/merge_adjacent_alignments.o : ${srcDir}/merge_adjacent_alignments.c ${libHeaders} abPOA
	${CC} ${CFLAGS} ${LDFLAGS} -o ${srcDir}/merge_adjacent_alignments.o -c ${srcDir}/merge_adjacent_alignments.c

${srcDir}/taf.o : ${srcDir}/taf.c ${libHeaders}
	${CC} ${CFLAGS} ${LDFLAGS} -o ${srcDir}/taf.o -c ${srcDir}/taf.c

${srcDir}/tai.o : ${srcDir}/tai.c ${libHeaders}
	${CC} ${CFLAGS} ${LDFLAGS} -o ${srcDir}/tai.o -c ${srcDir}/tai.c

${srcDir}/prefix_sort.o : ${srcDir}/prefix_sort.c ${libHeaders}
	${CC} ${CFLAGS} ${LDFLAGS} -o ${srcDir}/prefix_sort.o -c ${srcDir}/prefix_sort.c

${srcDir}/wiggle.o : ${srcDir}/wiggle.c ${libHeaders}
	${CC} ${CFLAGS} ${LDFLAGS} -o ${srcDir}/wiggle.o -c ${srcDir}/wiggle.c

${BINDIR}/stTafTests : ${libTests} ${LIBDIR}/libstTaf.a ${stTafDependencies}
	${CC} ${CFLAGS} ${LDFLAGS} -o ${BINDIR}/stTafTests ${libTests} ${LIBDIR}/libstTaf.a ${LDLIBS}

${BINDIR}/taffy : taf_norm.o taf_add_gap_bases.o taf_index.o taf_view.o taf_sort.o taf_stats.o taf_coverage.o taf_annotate.o taffy_main.o ${LIBDIR}/libstTaf.a ${libHeaders} ${stTafDependencies}
	${CXX} ${CPPFLAGS} ${CXXFLAGS} taf_norm.o taf_add_gap_bases.o taf_index.o taf_view.o taf_sort.o taf_stats.o taf_coverage.o taf_annotate.o taffy_main.o -o ${BINDIR}/taffy ${LIBDIR}/libstTaf.a ${LDLIBS}

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

taf_sort.o : taf_sort.c ${stTafDependencies} ${libHeaders}
	${CC} ${CFLAGS} ${CFLAGS} -o taf_sort.o -c taf_sort.c

taf_stats.o : taf_stats.c ${stTafDependencies} ${libHeaders}
	${CC} ${CFLAGS} ${CFLAGS} -o taf_stats.o -c taf_stats.c

taf_coverage.o : taf_coverage.cpp ${stTafDependencies} ${libHeaders}
	${CXX} ${CPPFLAGS} ${CXXFLAGS} -o taf_coverage.o -c taf_coverage.cpp

taf_annotate.o : taf_annotate.c ${stTafDependencies} ${libHeaders}
	${CC} ${CFLAGS} ${CFLAGS} -o taf_annotate.o -c taf_annotate.c

test : all ${BINDIR}/stTafTests
	${BINDIR}/stTafTests
	tests/tai/test_tai.py

python_test: all ${BINDIR}/stTafTests
	cd tests && python3 taffyTest.py

clean :
	cd ${sonLibRootDir} && ${MAKE} clean
	cd taffy/submodules/abPOA && ${MAKE} clean
	rm -rf *.o taffy/impl/*.o ${LIBDIR} ${BINDIR}

static :
	CFLAGS="$${CFLAGS} -static -march=nehalem" \
	CPPFLAGS="$${CXXFLAGS} -static -march=nehalem" \
	TAF_STATIC=1 \
	make all

python :
	python3 -m build
	python3 -m pip install .
	cd tests && python3 taffyTest.py && cd ..
