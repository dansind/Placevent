Oct 9, 2013
Upgrade to 1.4
1) Added Molecular Design Frontier 3D-RISM (MDF 3D-RISM) .h5 compatibility.
2) Use grid.py instead of internal pgrid module

Mar 23, 2013
Upgrade to 1.3
1) Steal population from pseudobulk (outside grid) instead of failing. 
Allows placement of many more atoms within smaller grid.

Mar 20, 2013
Upgrade to 1.2.2
1) Modified grid.py to read slightly differently formatted .dx files

Jan 28, 2013 
Upgrade to 1.2

1) Typo: Fixed minor error in equation for evacuation of final shell (will give slightly different results)
2) Adapted to setup.py-based installation 
    a) Added setup.py
    b) Rewrote tests to allow for either installed or unzipped placevent
    c) Uploaded to Pypi (try pip install placevent) 
       
