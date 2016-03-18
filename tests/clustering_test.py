#
"""
Tests for eqclustering module

- Mark C. Williams (2016) Nevada Seismological Lab

"""
import os
import logging
from eqclustering import BPTree

PWD = os.path.dirname(__file__)


def test_benchmark():
    """Benchmark reading/writing a file and output the time"""
    logging.captureWarnings(True)

    FINPUT = os.path.join(PWD, 'data', 'NSL_catalog_col.txt')
    FOUTPUT = "/tmp/eqclusterng-test-output.txt"
    
    import time
    t0 = time.time()
    
    t = BPTree.from_file(FINPUT) # Init tree with events from a file
    t.grow()                # Populate B-P tree with events
    t.prune()               # Prune given cutoff (calculate n if none)
    t.output2file(FOUTPUT)  # Output to file to match MATLAB perf

    t1 = time.time()
    print "Time: {0}s ".format(t1-t0),


