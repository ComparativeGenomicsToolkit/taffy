#!/usr/bin/env python3

#Copyright (C) 2006-2012 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt
import unittest
import bioioTest
import cigarsTest
from sonLib import treeTest
try:
    import networkx as NX
    networkx_installed = True
    from sonLib import nxtreeTest
    from sonLib import nxnewickTest
except ImportError:
    networkx_installed = False

from sonLib.bioio import system
from sonLib.bioio import parseSuiteTestOptions
from sonLib.bioio import getLogLevelString
from subprocess import check_call


def needsProgram(program):
    """Decorator: Run this test only if "program" is available."""
    def wrap(fn):
        try:
            check_call([program, "--version"])
        except:
            return unittest.skip(program + ' command is missing')(fn)
        else:
            return fn
    return wrap

def needsPackage(package):
    """Decorator: Run this test only if "pkg-config --exists package" says OK!"""
    def wrap(fn):
        try:
            check_call(["pkg-config", "--exists", package])
        except:
            return unittest.skip(package + ' pkg-config package is missing')(fn)
        else:
            return fn
    return wrap

class TestCase(unittest.TestCase):
    def testSonLibCTests(self):
        """Run most of the sonLib CuTests, fail if any of them fail."""
        system("sonLibTests %s" % getLogLevelString())

def allSuites():
    bioioSuite = unittest.makeSuite(bioioTest.TestCase, 'test')
    cigarsSuite = unittest.makeSuite(cigarsTest.TestCase, 'test')
    treeSuite = unittest.makeSuite(treeTest.TestCase, 'test')
    cuTestsSuite = unittest.makeSuite(TestCase, 'test')
    if not networkx_installed:
        allTests = unittest.TestSuite((bioioSuite, cigarsSuite, treeSuite, cuTestsSuite))
    else:
        nxtreeSuite = unittest.makeSuite(nxtreeTest.TestCase, 'test')
        nxnewickSuite = unittest.makeSuite(nxnewickTest.TestCase, 'test')
        allTests = unittest.TestSuite((bioioSuite, cigarsSuite, treeSuite, cuTestsSuite,
                                       nxtreeSuite, nxnewickSuite))
    return allTests

def main():
    parseSuiteTestOptions()

    suite = allSuites()
    runner = unittest.TextTestRunner(verbosity=2)
    i = runner.run(suite)
    return len(i.failures) + len(i.errors)

if __name__ == '__main__':
    import sys
    sys.exit(main())
