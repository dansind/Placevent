rm test/*pdb
echo "Running local placevent installation..."
python ../../placevent.py 1L2Y.O.1.dx 55.5 > test/O.pdb
echo "Printing difference between local test output and reference to file 'diffs'"
diff test/O.pdb ref/O.pdb > diffs
anydiffs=`wc -l diffs | awk '{print $1}'`;
if [ $anydiffs -gt 0 ] 
then
    echo "!! Differences found! Please check diffs file!"
    echo "Try downloading the latest version on: https://github.com/dansind/Placevent"
    echo "If this still fails, please contact me at sindhikara@gmail.com"
else 
    echo "No differences found! :)"
fi
