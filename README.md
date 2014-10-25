QMICAD
======

Quantum Mechanics Inspired Computer Aided Design (QMICAD) library is a 
collection of C++ classes and functions for simulation and design of 
nano-scaled devices using quantum mechanical tools such as Non-equilibrium 
Green Function (NEGF) formalism. This library provides a common framework 
for NEGF so that it can be used with any empirical tight binding, k.p model, 
extended Huckel method and density functional theory codes. This library also
includes a generic empirical tight binding model and a k.p model which can be 
extended for any material with known parameters.

Installation
-------------
Download and unzip source. Then configure using cmake:

`$ ./configure`

Make and install using:

`$ ./configure build`

`$ ./configure install ~/` 

To change installation directory, change `~/` to `/your/install/path`.

To create documentation do:

`$ ./configure build doc`

Copyright Notice
----------------

`Quantum Mechanics Inspired Computer Aided Design (QMICAD)
Copyright (C) 2014  K M Masum Habib <masum.habib@gmail.com

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.`
          
 
