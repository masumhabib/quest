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
Download and unzip source. Create a folder named 'build':
$ mkdir build
$ cd build
Then configure using cmake:
$ cmake ../
Make and install using:
$ make
$ make install DESTDIR=~/ 
To change installation directory, change DESTDIR option.
To create documentation do:
$ make doc

Copyright Notice
----------------

Quantum Mechanics Inspired Computer Aided Design (QMICAD)
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
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.
          
 
