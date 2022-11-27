#!/usr/bin/env python3

#Copyright (C) 2006-2012 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt

import unittest

import os
import sys
import random
import math
import io

from sonLib.bioio import getTempFile
from sonLib.bioio import getTempDirectory
from sonLib.bioio import TempFileTree
from sonLib.bioio import getRandomAlphaNumericString

from sonLib.bioio import fastaRead
from sonLib.bioio import fastaWrite
from sonLib.bioio import fastqRead
from sonLib.bioio import fastqWrite
from sonLib.bioio import getRandomSequence

from sonLib.bioio import pWMRead
from sonLib.bioio import pWMWrite

from sonLib.bioio import newickTreeParser
from sonLib.bioio import printBinaryTree

from sonLib.bioio import cigarRead
from sonLib.bioio import cigarWrite
from sonLib.bioio import PairwiseAlignment
from sonLib.bioio import getRandomPairwiseAlignment

from sonLib.bioio import TestStatus

from sonLib.bioio import system
from sonLib.bioio import logger

class TestCase(unittest.TestCase):

    def setUp(self):
        self.testNo = TestStatus.getTestSetup()
        unittest.TestCase.setUp(self)
        self.tempDir = getTempDirectory(os.getcwd())
        self.tempFiles = []

    def tearDown(self):
        unittest.TestCase.tearDown(self)
        system("rm -rf %s" % self.tempDir)
        for tempFile in self.tempFiles:
            os.remove(tempFile)

    #########################################################
    #########################################################
    #########################################################
    #temp file functions
    #########################################################
    #########################################################
    #########################################################

    def testTempFileTree(self):
        for test in range(100): #self.testNo):
            levels = random.choice(range(1, 4))
            fileNo = random.choice(range(1, 6))
            maxTempFiles = int(math.pow(fileNo, levels))

            print("Got %s levels, %s fileNo and %s maxTempFiles" % (levels, fileNo, maxTempFiles))

            tempFileTreeRootDir = os.path.join(self.tempDir, getRandomAlphaNumericString())
            tempFileTree = TempFileTree(tempFileTreeRootDir, fileNo, levels)

            tempFiles = []
            tempDirs = []
            #Check we can mac number of temp files.
            for i in range(maxTempFiles):
                if random.random() > 0.5:
                    tempFile = tempFileTree.getTempFile()
                    assert os.path.isfile(tempFile)
                    tempFiles.append(tempFile)
                else:
                    tempFile = tempFileTree.getTempDirectory()
                    assert os.path.isdir(tempFile)
                    tempDirs.append(tempFile)

            #Check assertion is created
            try:
                tempFileTree.getTempFile()
                assert False
            except RuntimeError:
                logger.debug("Got expected error message")

            #Now remove a few temp files
            while random.random() > 0.1 and len(tempFiles) > 0:
                tempFile = tempFiles.pop()
                assert os.path.isfile(tempFile)
                tempFileTree.destroyTempFile(tempFile)
                assert not os.path.isfile(tempFile)

            #Now remove a few temp dirs
            while random.random() > 0.1 and len(tempDirs) > 0:
                tempDir = tempDirs.pop()
                assert os.path.isdir(tempDir)
                tempFileTree.destroyTempDir(tempDir)
                assert not os.path.isdir(tempDir)

            #Check temp files is okay
            set(tempFileTree.listFiles()) == set(tempFiles + tempDirs)

            #Either remove all the temp files or just destroy the whole thing
            if random.random() > 0.5:
                #Remove all temp files and check thing is empty.
                for tempFile in tempFiles:
                    tempFileTree.destroyTempFile(tempFile)
                for tempDir in tempDirs:
                    tempFileTree.destroyTempDir(tempDir)
                os.remove(os.path.join(tempFileTreeRootDir, "lock"))
                os.rmdir(tempFileTreeRootDir)
            else:
                tempFileTree.destroyTempFiles()
                assert not os.path.isdir(tempFileTreeRootDir)

    #########################################################
    #########################################################
    #########################################################
    #fasta functions
    #########################################################
    #########################################################
    #########################################################

    def testFastaReadWrite(self):
        tempFile = getTempFile()
        self.tempFiles.append(tempFile)
        for test in range(0, self.testNo):
            fastaNumber = random.choice(range(10))
            l = [ getRandomSequence() for i in range(fastaNumber) ]
            fileHandle = open(tempFile, 'w')
            for name, seq in l:
                fastaWrite(fileHandle, name, seq)
            fileHandle.close()
            fileHandle = open(tempFile, 'r')
            l.reverse()
            outFh = io.StringIO()
            for i in fastaRead(fileHandle):
                assert i == l.pop()
                name, seq = i
                fastaWrite(outFh, name, seq)
            outFh.close()
            fileHandle.close()

    def testFastqReadWrite(self):
        tempFile = getTempFile()
        self.tempFiles.append(tempFile)
        for test in range(0, self.testNo):
            fastaNumber = random.choice(range(10))
            fastqs = [ (name, seq, [ random.randint(33, 126) for i in range(len(seq)) ]) for name, seq in [ getRandomSequence() for i in range(fastaNumber) ]]
            fH = open(tempFile, 'w')
            for name, seq, quals in fastqs:
                fastqWrite(fH, name, seq, quals)
            fH.close()
            fastqs.reverse()
            outFh = io.StringIO()
            for i in fastqRead(tempFile):
                assert i == fastqs.pop()
                name, seq, quals = i
                fastqWrite(outFh, name, seq, quals)
            outFh.close()

    def testFastaReadWriteC(self):
        """Tests consistency with C version of this function.
        """
        tempFile = getTempFile()
        self.tempFiles.append(tempFile)
        tempFile2 = getTempFile()
        self.tempFiles.append(tempFile2)
        for test in range(0, self.testNo):
            fastaNumber = random.choice(range(10))
            l = [ getRandomSequence() for i in range(fastaNumber) ]
            fileHandle = open(tempFile, 'w')
            for name, seq in l:
                fastaWrite(fileHandle, name, seq)
            fileHandle.close()

            command = "sonLib_fastaCTest %s %s" % (tempFile, tempFile2)

            print(command)

            system(command)

            fileHandle = open(tempFile2, 'r')
            l.reverse()
            outFh = io.StringIO()
            for i in fastaRead(fileHandle):
                name, seq = i
                assert i == l.pop()
                fastaWrite(outFh, name, seq)
            outFh.close()
            fileHandle.close()

    #########################################################
    #########################################################
    #########################################################
    #newick tree functions
    #########################################################
    #########################################################
    #########################################################

    def testNewickTreeParser(self):
        if self.testNo > 0:
            d = '((human,baboon),chimp);'
            e = newickTreeParser(d)
            f = printBinaryTree(e, False)
            print(d, f)
            assert d == f

    def testNewickTreeParser_UnaryNodes(self):
        #tests with unary nodes
        for test in range(0, self.testNo):
            tree = getRandomTreeString()
            logger.debug("tree to try\t", tree)
            tree2 = newickTreeParser(tree, reportUnaryNodes=True)
            tree3 = printBinaryTree(tree2, True)
            logger.debug("tree found\t", tree3)
            assert tree == tree3

    #########################################################
    #########################################################
    #########################################################
    #functions for position weight matrices
    #########################################################
    #########################################################
    #########################################################

    def testPWMParser(self):
        tempFile = getTempFile()
        self.tempFiles.append(tempFile)
        for test in range(0, self.testNo):
            pWM = getRandomPWM()

            fileHandle = open(tempFile, 'w')
            pWMWrite(fileHandle, pWM)
            fileHandle.close()

            fileHandle = open(tempFile, 'r')
            pWM2 = pWMRead(fileHandle)
            fileHandle.close()

            for i in range(0, len(pWM)):
                pWM[i] == pWM2[i]

    #########################################################
    #########################################################
    #########################################################
    #Cigar functions
    #########################################################
    #########################################################
    #########################################################

    def testCigarReadWrite(self):
        tempFile = getTempFile()
        self.tempFiles.append(tempFile)
        for test in range(0, self.testNo):
            cigarNumber = random.choice(range(10))
            l = [ getRandomPairwiseAlignment() for i in range(cigarNumber) ]
            fileHandle = open(tempFile, 'w')
            for cigar in l:
                cigarWrite(fileHandle, cigar)
            fileHandle.close()
            fileHandle = open(tempFile, 'r')
            l.reverse()
            for cigar in cigarRead(fileHandle):
                cigarWrite(sys.stdout, l[-1])
                cigarWrite(sys.stdout, cigar)
                assert cigar == l.pop()
            assert len(l) == 0
            fileHandle.close()

def getRandomTreeString():
    def iDFn():
        return random.choice([ "one", "1", "", "he44"])
    def dFn():
        #if random.random() > 0.5:
        return ":%.6f" % random.random()
        #return ''
    def fn3():
        if random.random() > 0.5:
            if random.random() > 0.5:
                #is unary
                return '(' + fn3() + ')' + iDFn() + dFn()
            else:
                return '(' + fn3() +  ","  + fn3() + ')' + iDFn() + dFn()
        else:
            return iDFn() + dFn()
    return fn3() + ';'

def getRandomPWM(length=-1):
    if length == -1:
        length = 1 + int(random.random()*10)
    def fn():
        l = [ random.random()*100 for i in range(0, 4) ]
        i = sum(l)
        return [ j/i for j in l ]
    return [ fn() for i in range(0, length) ]

if __name__ == '__main__':
    unittest.main()
