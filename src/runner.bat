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

GOTO Jmol

GOTO Finished



:MultiDR

echo Testing Multiple DR Files

python TBModeller.py testdata\T0644.dr testdata\T0645.dr testdata\T0648.dr

pause

GOTO Finished



:Jmol

echo Visualizing Model

REM Start Jmol

java -Xmx512m -jar "C:\Jmol\Jmol.jar" targets\T0644_chain.pdb

GOTO Finished



REM Done!

:Finished

echo Finished!
