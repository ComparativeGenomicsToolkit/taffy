include ../../include.mk
BINDIR = ../../bin
LIBDIR = ../../lib

all : ${LIBDIR}/cuTest.a

clean : 
	rm -f ${LIBDIR}/cuTest.a ${LIBDIR}/CuTest.h *.o

${LIBDIR}/cuTest.a : CuTest.c CuTest.h
	@mkdir -p $(dir $@)
	${CC} ${CFLAGS} -c CuTest.c
	${AR} rc cuTest.a *.o
	${RANLIB} cuTest.a
	mv cuTest.a ${LIBDIR}/
	cp CuTest.h ${LIBDIR}/
