eqclustering
============
Python library implementation of earthquake clustering algorithm of Zaliapin et al.

This is a python implmentation of the clustering algortihm, based on the original MATLAB implementation of Zaliapin. The only hard python lib dependency is numpy. It is recommended to use a numpy compiled against ATLAS, BLAS, LAPACK, etc. for best performance. 

WORK IN PROGRESS!!
------------------
This software is in flux. It currently produces the same results as the MATLAB code used for the original publication...

...But the I/O methods, tests and test data, and other features may change rapidly until a 1.0 release.

* Documentation - here and improved docstrings
* I/O - make more agnostic
* numpy - better integration, allow for using numpy.load, use 2D array with attr, etc

Installation
------------
Recommended install is via pip

Dependencies
------------
numpy - Use a numpy built against ATLAS, BLAS, LAPACK for optimal performance

Citation
--------
The authors request any citation, especially published results using this software, include the author and publication of the original algorithm (Zaliapin):

 Zaliapin et al. (2008) "Clustering Analysis of Seismicity and Aftershock Identification" Phys. Rev. Lett. 101, 018501 
 doi:10.1103/PhysRevLett.101.018501

 University of Nevada, Reno (1971): Nevada Seismic Network. International Federation of Digital Seismograph Networks. Other/Seismic Network.
 doi:10.7914/SN/NN

Authors
-------
Mark Williams (Nevada Seismological Laboratory)

Copyright
---------
Copyright 2016 University of Nevada, Reno

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

   http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

