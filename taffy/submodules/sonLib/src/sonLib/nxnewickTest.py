#!/usr/bin/env python3

#Copyright (C) 011 by Glenn Hickey
#
#Released under the MIT license, see LICENSE.txt

import unittest
import random
import re
from string import whitespace as ws

from sonLib.nxnewick import NXNewick
from sonLib.bioio import TestStatus
from sonLib.bioio import logger

class TestCase(unittest.TestCase):
    
    def setUp(self):
        self.testNo = TestStatus.getTestSetup()
        unittest.TestCase.setUp(self)
    
    def tearDown(self):
        unittest.TestCase.tearDown(self)
    
    def testNewickIO(self):
        # feslenstein's own... (http://evolution.genetics.washington.edu/phylip/newicktree.html)
        tree1 = '((raccoon:19.19959,bear:6.80041):0.846,((sea_lion:11.997, seal:12.003):7.52973,((monkey:100.8593,cat:47.14069):20.59201, weasel:18.87953):2.0946):3.87382,dog:25.46154);'
        tree2 = '(Bovine:0.69395,(Gibbon:0.36079,(Orang:0.33636,(Gorilla:0.17147,(Chimp:0.19268, Human:0.11927):0.08386):0.06124):0.15057):0.54939,Mouse:1.2146):0.1;'
        tree3 = '(Bovine:0.69395,(Hylobates:0.36079,(Pongo:0.33636,(G._Gorilla:0.17147, (P._paniscus:0.19268,H._sapiens:0.11927):0.08386):0.06124):0.15057):0.54939, Rodent:1.2146);'
        tree4 = 'A;'
        tree5 = '((A,B):0.0,(C,D));'
        tree6 = '(Alpha,Beta,Gamma,Delta,,Epsilon,,,);'
        
        trees = [tree1, tree2, tree3, tree4, tree5, tree6]        
        newickParser = NXNewick()
        
        # Parse newicks, adding implied roots
        for tree in trees:
            newickParser.parseString(tree, addImpliedRoots=True)
            answer = self.__cleanTree(tree)
            outputString = newickParser.writeString()
            logger.debug(" ***************** ")
            logger.debug(outputString)
            logger.debug(answer)
            assert outputString == answer
    
        # Parse newicks, not adding implied roots
        for tree in trees:
            newickParser.parseString(tree, addImpliedRoots=False)
            outputString = newickParser.writeString()
            answer = re.sub(r':[.0-9]+?;', ';', tree)
            answer = re.sub(r'\s+', '', answer)
            logger.debug(" ***************** ")
            logger.debug(outputString)
            logger.debug(answer)
            assert outputString == answer
        

    # remove whitespace
    # change trees of form (blabla):123; to ((blabla):123);
    def __cleanTree(self, tree):
        result = ""
        for c in tree:
            if c not in ws: 
                result += c
        if result[len(result) -2 ] != ')' and result[0] == '(':
                result = "(%s);" % result[:len(result) - 1]
        return result
    
if __name__ == '__main__':
    unittest.main()
