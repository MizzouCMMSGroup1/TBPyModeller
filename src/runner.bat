@echo OFF

REM TBModeller Test Runner

echo Testing Single Sequence

python TBModeller.py T0644 MKFLKFSLLTAVLLSVVFAFSSCGDDDDTGYLPPSQAIQDALKKLYPNATAIKWEQKGVYYVADCQADGREKEVWFDANANWLMTETELNSINNLPPAVLTAFMESSYNNWVVDDVVILEYPNEPSTEFVVTVEQGKKVDLYFSEGGGLLHEKDVTNGDDTHWPRV

pause

GOTO Finished

echo Testing DR File

python TBModeller.py testdata\T0644.dr

pause

echo Testing FASTA File

python TBModeller.py testdata\T0644.fasta

pause

echo Testing Multiple DR Files

python TBModeller.py testdata\T0644.dr testdata\T0645.dr testdata\T0648.dr

:Finished

echo Finished!
