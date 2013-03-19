echo "Running placevent..."
if [ -x "`which placevent.py`" ] ; then
    exampledir=`python -c "import placevent; import os; print(os.path.join(os.path.dirname(placevent.__file__),'tests','1L2Y'))"`
    placevent.py $exampledir/1L2Y.O.1.dx 55.5 | grep "ATOM" > O_test.pdb
    diff O_test.pdb $exampledir/ref/O.pdb > diffs
else
    #Assuming this is being run in the test directory
    ../../placevent.py 1L2Y.O.1.dx 55.5 | grep "ATOM" > O_test.pdb
    diff O_test.pdb ref/O.pdb > diffs
fi

echo "Printing difference between local test output and reference to file 'diffs'"
anydiffs=`wc -l diffs | awk '{print $1}'`;
if [ $anydiffs -gt 0 ] 
then
    echo "!! Differences found! Please check diffs file!"
    echo "Try downloading the latest version on: https://github.com/dansind/Placevent"
    echo "If this still fails, please contact me at sindhikara@gmail.com"
else 
    echo "No differences found! :)"
fi

