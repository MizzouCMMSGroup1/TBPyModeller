#/bin/bash



echo 'Running TBModeller.py';
echo 'Testing Single Sequence';

python3 TBModeller.py T0644 MKFLKFSLLTAVLLSVVFAFSSCGDDDDTGYLPPSQAIQDALKKLYPNATAIKWEQKGVYYVADCQADGREKEVWFDANANWLMTETELNSINNLPPAVLTAFMESSYNNWVVDDVVILEYPNEPSTEFVVTVEQGKKVDLYFSEGGGLLHEKDVTNGDDTHWPRV;

echo 'Testing DR File';

python3 TBModeller.py 'testdata/T0644.dr';

echo "Testing FASTA File"

python3 TBModeller.py "testdata/T0644.fasta";

echo "Testing Multiple DR Files"

python3 TBModeller.py "testdata/T0644.dr" "testdata/T0645.dr" "testdata/T0648.dr"

echo 'Finished!'
