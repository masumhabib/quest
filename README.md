QUEST
======

QUantum mechanics Enabled Simulation Toolset (QUEST) is a 
collection of C++ classes and functions for simulation and design of 
nano-scaled devices using quantum mechanical tools such as Non-equilibrium 
Green Function (NEGF) formalism. This library provides a common framework 
for NEGF so that it can be used with any empirical tight binding, k.p model, 
extended Huckel method and density functional theory codes. This library also
includes a generic empirical tight binding model and a k.p model which can be 
extended for any material with known parameters.

Unlike most of the quantum transport tools available out there, 
QUEST is a *not* a predefined simulator for a predefined problem in 
a predefined device structure. Rather, it is a flexible toolbox 
that allows users to write their own simulators best-suited for their 
own problems. 

Deep down, QUEST is a C++ library that uses BLAS and LAPACK for computation.
Yet, QUEST provides both C++ and Python interfaces so that one can choose between 
the elegance of C++ and ease of Python.

Dependencies
-------------

The following libraries are needed to build QUEST:

* CMake >= 2.8
* Boost >= 1.55 (Boost.MPI and Boost.Python)
* Armadillo >= 4.x
* BLAS+LAPACK (Intel MKL)
* Python >= 2.7
* Doxygen >= 1.8

Installation
-------------
Download and unzip source. Then configure using cmake:   
`$ ./build config`

Make and install using:   
`$ ./build`   
`$ ./build install ~/`

To change installation directory, change `~/` to `/your/install/path`.

To create documentation do:   
`$ ./build doc`


Documentation
--------------

Please see the documentation produced by doxygen. 
We also have a [wiki](https://github.com/masumhabib/quest/wiki/Home).

Copyright Notice
----------------
QUantum mechanics Enabled Simulation Toolset (QUEST)
Copyright (C) 2013-2016  by the authors listed in the NOTICE file. 

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
