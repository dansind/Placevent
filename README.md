Placevent
=========

Placevent - 3D-RISM-based solvent and ion placement software created by Daniel Sindhikara, sindhikara@gmail.com
Placevent
This program is designed to automatically place explicit solvent atoms/ions based
on 3D-RISM data. Placevent is compatible with AMBER .dx files and MDF .h5 files.

The details of the algorithm are described in: 
Placevent: An algorithm for prediction of explicit solvent atom distribution -- Awpplication to HIV-1 protease and F-ATP synthase 
Daniel J. Sindhikara, Norio Yoshida, Fumio Hirata
http://dansindhikara.com/Software/Entries/2012/6/22_Placevent_New.html

If you have Numpy, this package should work out of the box. Please contact me if you have problems.
This package requires Numpy.
This package requires grid.py.

EXAMPLE:
Please see the tests/1L2Y/ for examples/tests
In that directory, run test.sh to see if your local installation matches the reference output.

TUTORIAL:
See my website: http://dansindhikara.com/Tutorials/Entries/2012/1/1_Using_3D-RISM_and_PLACEVENT.html
for a tutorial.


g(r)_0 is printed in the occupation column
g(r)_i is printed in the beta column



    Copyright (C) 2012 Daniel J. Sindhikara

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
