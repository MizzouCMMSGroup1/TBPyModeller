#/bin/bash



echo 'Running TBModeller.py';
echo 'Testing Single Sequence';

python TBModeller.py T0644 MKFLKFSLLTAVLLSVVFAFSSCGDDDDTGYLPPSQAIQDALKKLYPNATAIKWEQKGVYYVADCQADGREKEVWFDANANWLMTETELNSINNLPPAVLTAFMESSYNNWVVDDVVILEYPNEPSTEFVVTVEQGKKVDLYFSEGGGLLHEKDVTNGDDTHWPRV;

echo 'Testing DR File';

python TBModeller.py 'testdata/T0644.dr';

echo "Testing FASTA File"

python TBModeller.py "testdata/T0644.fasta";

echo "Testing Multiple DR Files"

python TBModeller.py "testdata/T0644.dr" "testdata/T0645.dr" "testdata/T0648.dr"

echo 'Finished!'
