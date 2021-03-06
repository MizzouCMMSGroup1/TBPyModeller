TBPyModeller
============

Template Based Modelling of Protein Structures in Python

#####Project coordinator:  

| Name          | Handle
|---------------|--------------------------------------
| Xiaokai Qian  | [qxk2310](https://github.com/qxk2310)

#####Team Members:  

| Name          | Handle                                     
|---------------|------------------------------------------
| Caiwei Wang   | [cwwmoment](https://github.com/cwwmoment)
| Haipei Fan    | [happyfhp](https://github.com/happyfhp)
| Puneet Gaddam | [pg7b2](https://github.com/pg7b2)
| Sean Lander   | [slander36](https://github.com/slander36)
| Brett Koonce  | [asparagui](https://github.com/asparagui)

Requirements:
=============
#####Required:  
* [Python 3.3 32 bit](http://www.python.org/download/)  
* [Numpy 1.8 for Python 3.3 32 bit](http://www.numpy.org/)  
* [Biopython 1.63 for Python 3.3 32 bit](http://biopython.org/wiki/Download)  
* [Scwrl4](http://dunbrack.fccc.edu/scwrl4/)  

#####Optional:  
* [Jmol](http://jmol.sourceforge.net/download/)  

To Run:
=======
```
# Using ID and Sequence  
$ python TBModeller.py ID Sequence [ID Sequence...]  

# Using Fasta Files
$ python TBModeller.py file1.fasta [file2.fasta...]  

# Using DR Files (CASP)
$ python TBModeller.py file1.dr [file2.dr...]  

# Visualization with Jmol on Windows
C:\> set JMOL_HOME="{Jmol Directory Here}"
# Example: C:\> set JMOL_HOME="C:\Jmol"   
C:\> java -Xmx512m -jar %JMOL_HOME%\Jmol.jar targets\{targetid}.pdb

# Visualization with Jmol on Linux/Mac
$ export JMOL_HOME="{Jmol Directory Here}"
# Example: $ export JMOL_HOME="/home/user/Jmol"
$ java -Xmx512m -jar "$JMOL_HOME/Jmol.jar" targets/{targetid}.pdb
```
