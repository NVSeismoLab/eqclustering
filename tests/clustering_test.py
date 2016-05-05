# -*- coding: utf-8 -*-
#
# Copyright 2016 University of Nevada, Reno
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# Author: Mark Williams, Nevada Seismological Laboratory, UNR
#
"""
Tests for eqclustering module

- Mark C. Williams (2016) Nevada Seismological Lab

"""
import os
import logging
from eqclustering import BPTree

PWD = os.path.dirname(__file__)


def test_bptree():
    """Test creating and pruning tree from input data"""
    logging.captureWarnings(True)

    FINPUT = os.path.join(PWD, 'data', 'NSL_catalog_col.txt')
    
    t = BPTree.from_file(FINPUT)
    t.grow()
    assert len(t.P) == 283  # check number of events in file
    
    t.prune()
    assert len(t.Pfin) == 283    


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
    print " bench: {0}eq/{1:.6f}s ".format(len(t.P), t1-t0),


