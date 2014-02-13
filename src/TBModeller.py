"""

TBModeller - A simple Template-Based Modelling program for Protein Sequences

Copyright (c) Sean Lander 2014 under GPL v3

"""

from Protein import Protein
from Bio import SeqIO
from Bio.Seq import Seq
import sys
import re

targets = []

def main():
	for target in targets:
		print(target)

if __name__ == "__main__":
	# RUN ME
	# Check if sequence or fasta
	if len(sys.argv) != 2:
		print("Error: Expected 2 arguments, received {0}".format(len(sys.argv)))
		sys.exit()
	input = sys.argv[1]
	if ".fasta" in input:
		# Fasta file
		for seq_record in SeqIO.parse(input, "fasta"):
			targets.append(Protein(seq_record.seq))
	elif ".dr" in input:
		# DR file
		print("Opening DR file")
		f = open(input,'r')
		sequence = []
		for line in f:
			# Check format, use residue lines only
			if not re.search("\w\s\w\s\d+",line):
				continue
			sequence.append(line.split('\t')[0])
		sequence = ''.join(sequence)
		targets.append(Protein(Seq(sequence)))
	else:
		# Raw sequence
		targets.append(Protein(Seq(input)))

	main()