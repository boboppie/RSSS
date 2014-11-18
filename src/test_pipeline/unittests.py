#!/usr/bin/env python

import sys; sys.path.append("core")
from steps import *

if __name__ == '__main__':
    #unittest.main(exit=True, verbosity=2)
    print "Start unit testing ..."
    print "================ [SimulatedReadsTests] =================="
    SimulatedReadsTests()

    print "Unit tests finished"
