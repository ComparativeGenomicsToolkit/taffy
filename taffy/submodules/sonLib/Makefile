include include.mk
BINDIR = ./bin

.PHONY: all clean cP cP.clean externalToolsP.clean test

all : cP ${BINDIR}/sonLib_daemonize.py

clean : cP.clean externalToolsP.clean
	rm -f ${BINDIR}/sonLib_daemonize.py

cP: externalToolsP
	cd C && ${MAKE} all

externalToolsP:
	cd externalTools && ${MAKE} all

cP.clean :
	cd C && ${MAKE} clean

externalToolsP.clean :
	cd externalTools && ${MAKE} clean

test : all
	PYTHONPATH=src:. PATH=$$(pwd)/bin:$$PATH ${PYTHON} allTests.py --testLength=SHORT --logLevel=WARN

${BINDIR}/sonLib_daemonize.py : sonLib_daemonize.py cP
	cp sonLib_daemonize.py ${BINDIR}/sonLib_daemonize.py
	chmod +x ${BINDIR}/sonLib_daemonize.py
