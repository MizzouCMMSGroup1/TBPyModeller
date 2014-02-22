"""

TBModeller - A simple Template-Based Modelling program for Protein Sequences

Copyright (c) Sean Lander 2014 under GPL v2

"""

from Protein import Protein

from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import sys
import re
import os
import shutil

targets = []

debug = True

skipBlanks = False

def main():
	for target in targets:
		print("Target:",target)
		print()
		print(target.getTemplates())
		print(target.getPDBs())
		target.alignPDB(skipBlanks = skipBlanks)

		# If we are skipping blanks we can use Scwrl4
		# Otherwise, fake the operation for now
		if skipBlanks:
			os.system("Scwrl4 -i targets/%s.pdb -o targets/%s_chain.pdb" % (target,target))
		else:
			shutil.copyfile('targets/%s.pdb' % (target), 'targets/%s_chain.pdb' % (target))

if __name__ == "__main__":
	# RUN ME
	# Check if sequence or fasta
	if len(sys.argv) < 2:
		print("Error: Expected 2+ arguments, received {0}".format(len(sys.argv)))
		sys.exit()
	i = 1
	while i < len(sys.argv):

		input = sys.argv[i]
		if ".fasta" in input:
			# Fasta file
			for seq_record in SeqIO.parse(input, "fasta"):
				seq_record.seq.alphabet = IUPAC.protein
				targets.append(Protein(seq_record.id,seq_record.seq,debug=debug))
		elif ".dr" in input:
			# DR file
			f = open(input,'r')
			id = 0
			sequence = []
			for line in f:
				# Check format, use residue lines only
				if re.match("TARGET",line):
					id = line.strip().split('\t')[1].strip()
				if not re.search("\w\s\w\s\d+",line):
					continue
				sequence.append(line.split('\t')[0])
			sequence = ''.join(sequence)
			targets.append(Protein(id,Seq(sequence,IUPAC.protein),debug=debug))
		else:
			# Raw sequence
			id = input
			i = i+1
			input = sys.argv[i]
			targets.append(Protein(id,Seq(input,IUPAC.protein),debug=debug))
		i = i+1

	main()