rm test/*pdb
python ../../placevent.py 1L2Y.O.1.dx 55.5 > test/O.pdb
diff test/O.pdb ref/O.pdb
