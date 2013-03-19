===========
Placevent
===========
The details of the algorithm are described in: 

Placevent: An algorithm for prediction of explicit solvent atom distribution --
Application to HIV-1 protease and F-ATP synthase 
Daniel J. Sindhikara, Norio Yoshida, Fumio Hirata


Utilizes 3D-RISM output to place solvent or ions in highest likelihood locations.
Matches well with experiment in test cases.


Typical usage often looks like this::

    placevent.py mydxfile.dx 55.5 > O.pdb


