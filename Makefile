rootPath = .
include ${rootPath}/include.mk

srcDir = taffy/impl
oneCodeDir = taffy/submodules/ONEcode
bigWigDir = taffy/submodules/libBigWig
blockVizDir = taffy/blockViz
libHeaders = taffy/inc/*.h ${blockVizDir}/inc/*.h
libTests = tests/*.c

all: all_libs all_progs
all_libs: ${LIBDIR}/libstTaf.a

all_progs: all_libs ${BINDIR}/taffy ${BINDIR}/taffyBlockVizTest

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

# Vendored libBigWig (local-file reader only; -DNOCURL drops the libcurl dep).
# Built like abPOA into ${LIBDIR}/libBigWig.a and linked via -lBigWig.
${LIBDIR}/libBigWig.a :
	mkdir -p ${LIBDIR} ${INCLDIR}
	cd ${bigWigDir} && ${CC} -O3 -g -fPIC -DNOCURL -c bwRead.c bwStats.c bwValues.c bwWrite.c io.c
	${AR} rc ${LIBDIR}/libBigWig.a ${bigWigDir}/bwRead.o ${bigWigDir}/bwStats.o ${bigWigDir}/bwValues.o ${bigWigDir}/bwWrite.o ${bigWigDir}/io.o
	ln -f ${bigWigDir}/bigWig.h ${bigWigDir}/bwCommon.h ${bigWigDir}/bwValues.h ${bigWigDir}/bigWigIO.h ${INCLDIR}

${sonLibDir}/sonLib.a : sonLib

${sonLibDir}/cuTest.a : sonLib

stTafDependencies = ${sonLibDir}/sonLib.a ${sonLibDir}/cuTest.a ${LIBDIR}/libabpoa.a ${LIBDIR}/libBigWig.a

${oneCodeDir}/ONElib.o : ${oneCodeDir}/ONElib.c ${oneCodeDir}/ONElib.h
	${CC} ${CFLAGS} ${LDFLAGS} -fPIC -o ${oneCodeDir}/ONElib.o -c ${oneCodeDir}/ONElib.c

${LIBDIR}/libstTaf.a : ${libTests} ${libHeaders} ${srcDir}/alignment_block.o ${srcDir}/line_iterator.o ${srcDir}/maf.o ${srcDir}/paf.o ${srcDir}/ond.o ${srcDir}/taf.o ${srcDir}/add_gap_bases.o ${srcDir}/merge_adjacent_alignments.o ${srcDir}/prefix_sort.o ${srcDir}/wiggle.o ${srcDir}/tai.o ${srcDir}/tui.o ${srcDir}/chain.o ${srcDir}/view_chain_dup_filter.o ${srcDir}/block_reader.o ${srcDir}/remote_io.o ${srcDir}/gerp.o ${srcDir}/gerp_stats.o ${blockVizDir}/impl/taffyBlockViz.o ${oneCodeDir}/ONElib.o ${libHeaders} ${stTafDependencies}
	${AR} rc libstTaf.a ${srcDir}/alignment_block.o ${srcDir}/line_iterator.o ${srcDir}/maf.o ${srcDir}/paf.o ${srcDir}/ond.o ${srcDir}/taf.o ${srcDir}/add_gap_bases.o ${srcDir}/merge_adjacent_alignments.o ${srcDir}/prefix_sort.o ${srcDir}/wiggle.o ${srcDir}/tai.o ${srcDir}/tui.o ${srcDir}/chain.o ${srcDir}/view_chain_dup_filter.o ${srcDir}/block_reader.o ${srcDir}/remote_io.o ${srcDir}/gerp.o ${srcDir}/gerp_stats.o ${blockVizDir}/impl/taffyBlockViz.o ${oneCodeDir}/ONElib.o
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

${blockVizDir}/impl/taffyBlockViz.o : ${blockVizDir}/impl/taffyBlockViz.cpp ${blockVizDir}/inc/taffyBlockViz.h ${libHeaders}
	${CXX} ${CPPFLAGS} ${CXXFLAGS} -I${blockVizDir}/inc -o ${blockVizDir}/impl/taffyBlockViz.o -c ${blockVizDir}/impl/taffyBlockViz.cpp

${srcDir}/merge_adjacent_alignments.o : ${srcDir}/merge_adjacent_alignments.c ${libHeaders} abPOA
	${CC} ${CFLAGS} ${LDFLAGS} -o ${srcDir}/merge_adjacent_alignments.o -c ${srcDir}/merge_adjacent_alignments.c

${srcDir}/taf.o : ${srcDir}/taf.c ${libHeaders}
	${CC} ${CFLAGS} ${LDFLAGS} -o ${srcDir}/taf.o -c ${srcDir}/taf.c

${srcDir}/tai.o : ${srcDir}/tai.c ${libHeaders}
	${CC} ${CFLAGS} ${LDFLAGS} -o ${srcDir}/tai.o -c ${srcDir}/tai.c

${srcDir}/tui.o : ${srcDir}/tui.c ${libHeaders} ${oneCodeDir}/ONElib.o
	${CC} ${CFLAGS} ${LDFLAGS} -o ${srcDir}/tui.o -c ${srcDir}/tui.c

${srcDir}/chain.o : ${srcDir}/chain.c ${libHeaders}
	${CC} ${CFLAGS} ${LDFLAGS} -o ${srcDir}/chain.o -c ${srcDir}/chain.c

${srcDir}/view_chain_dup_filter.o : ${srcDir}/view_chain_dup_filter.c ${libHeaders}
	${CC} ${CFLAGS} ${LDFLAGS} -o ${srcDir}/view_chain_dup_filter.o -c ${srcDir}/view_chain_dup_filter.c

${srcDir}/prefix_sort.o : ${srcDir}/prefix_sort.c ${libHeaders}
	${CC} ${CFLAGS} ${LDFLAGS} -o ${srcDir}/prefix_sort.o -c ${srcDir}/prefix_sort.c

${srcDir}/block_reader.o : ${srcDir}/block_reader.c ${libHeaders}
	${CC} ${CFLAGS} ${LDFLAGS} -o ${srcDir}/block_reader.o -c ${srcDir}/block_reader.c

${srcDir}/wiggle.o : ${srcDir}/wiggle.c ${libHeaders}
	${CC} ${CFLAGS} ${LDFLAGS} -o ${srcDir}/wiggle.o -c ${srcDir}/wiggle.c

${srcDir}/remote_io.o : ${srcDir}/remote_io.c ${libHeaders}
	${CC} ${CFLAGS} ${LDFLAGS} -o ${srcDir}/remote_io.o -c ${srcDir}/remote_io.c

${BINDIR}/stTafTests : ${libTests} ${LIBDIR}/libstTaf.a ${stTafDependencies}
	${CC} ${CFLAGS} ${LDFLAGS} -I${blockVizDir}/inc -o ${BINDIR}/stTafTests ${libTests} ${LIBDIR}/libstTaf.a ${LDLIBS}

${BINDIR}/taffy : taf_norm.o taf_add_gap_bases.o taf_index.o taf_view.o taf_sort.o taf_stats.o taf_coverage.o taf_summary.o taf_annotate.o taf_lift.o taf_chain.o taf_gerp.o taf_gerp_stats.o taffy_main.o ${LIBDIR}/libstTaf.a ${libHeaders} ${stTafDependencies}
	${CXX} ${CPPFLAGS} ${CXXFLAGS} taf_norm.o taf_add_gap_bases.o taf_index.o taf_view.o taf_sort.o taf_stats.o taf_coverage.o taf_summary.o taf_annotate.o taf_lift.o taf_chain.o taf_gerp.o taf_gerp_stats.o taffy_main.o -o ${BINDIR}/taffy ${LIBDIR}/libstTaf.a ${LDLIBS}

${BINDIR}/taffyBlockVizTest : ${blockVizDir}/tests/taffyBlockVizTest.cpp ${LIBDIR}/libstTaf.a ${blockVizDir}/inc/taffyBlockViz.h ${stTafDependencies}
	${CXX} ${CPPFLAGS} ${CXXFLAGS} -I${blockVizDir}/inc ${blockVizDir}/tests/taffyBlockVizTest.cpp -o ${BINDIR}/taffyBlockVizTest ${LIBDIR}/libstTaf.a ${LDLIBS}

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

taf_summary.o : taf_summary.cpp ${stTafDependencies} ${libHeaders}
	${CXX} ${CPPFLAGS} ${CXXFLAGS} -I${blockVizDir}/inc -o taf_summary.o -c taf_summary.cpp

taf_annotate.o : taf_annotate.c ${stTafDependencies} ${libHeaders}
	${CC} ${CFLAGS} ${CFLAGS} -o taf_annotate.o -c taf_annotate.c

taf_lift.o : taf_lift.c ${stTafDependencies} ${libHeaders}
	${CC} ${CFLAGS} ${CFLAGS} -o taf_lift.o -c taf_lift.c

taf_chain.o : taf_chain.c ${stTafDependencies} ${libHeaders}
	${CC} ${CFLAGS} ${CFLAGS} -o taf_chain.o -c taf_chain.c

taf_gerp.o : taf_gerp.c ${stTafDependencies} ${libHeaders}
	${CC} ${CFLAGS} ${CFLAGS} -o taf_gerp.o -c taf_gerp.c

taf_gerp_stats.o : taf_gerp_stats.c ${stTafDependencies} ${libHeaders}
	${CC} ${CFLAGS} ${CFLAGS} -o taf_gerp_stats.o -c taf_gerp_stats.c

test : all ${BINDIR}/stTafTests
	${BINDIR}/stTafTests
	pytest tests/tai tests/tui

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
