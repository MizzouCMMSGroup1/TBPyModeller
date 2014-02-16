A quick demo of homology modelling using modeller
===

mostly cribbed from: http://salilab.org/modeller/tutorial/basic.html

short version:

download http://salilab.org/modeller/downloads/pdball.pir.gz to get the master sequence database

all command below via modeller (in other words run:  mod9.13 command.py)

build_profile.py
==

builds a local binary pdb database from the pdball file you downloaded, then

searches (BLOSUM) the master pdb database to get a list of similar sequence templates

compare.py
==
uses above profile to match the alignments of the found sequences (match the alignments)

then, I just picked a good template (not going to do combinations yet): 3U1W and made that the only template we're using

align2d_mod.py
==

uses the alignment from prior to build a T0644 PDB file from our input template (3U1W) to a new alignment file

build_model.py
==

generates a collection of models from the new alignment file

