@echo OFF

REM TBModeller Test Runner

echo Running TBModeller.py

GOTO SingleFasta



:SingleSeq

echo Testing Single Sequence

python TBModeller.py T0644 MKFLKFSLLTAVLLSVVFAFSSCGDDDDTGYLPPSQAIQDALKKLYPNATAIKWEQKGVYYVADCQADGREKEVWFDANANWLMTETELNSINNLPPAVLTAFMESSYNNWVVDDVVILEYPNEPSTEFVVTVEQGKKVDLYFSEGGGLLHEKDVTNGDDTHWPRV

pause

GOTO Finished



:SingleDR

echo Testing DR File

python TBModeller.py testdata\T0644.dr

pause

GOTO Finished



:SingleFasta

echo Testing FASTA File

python TBModeller.py testdata\T0644.fasta

pause

GOTO Finished



:MultiDR

echo Testing Multiple DR Files

python TBModeller.py testdata\T0644.dr testdata\T0645.dr testdata\T0648.dr

pause

GOTO Finished



REM Done!

:Finished

echo Finished!
