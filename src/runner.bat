@echo OFF

REM TBModeller Test Runner

echo Testing Single Sequence

python TBModeller.py MKFLKFSLLTAVLLSVVFAFSSCGDDDDTGYLPPSQAIQDALKKLYPNATAIKWEQKGVYYVADCQADGREKEVWFDANANWLMTETELNSINNLPPAVLTAFMESSYNNWVVDDVVILEYPNEPSTEFVVTVEQGKKVDLYFSEGGGLLHEKDVTNGDDTHWPRV

pause

echo Testing DR File

python TBModeller.py testdata\T0644.dr

pause

echo Testing FASTA File

python TBModeller.py testdata\T0644.fasta

pause

echo Finished!