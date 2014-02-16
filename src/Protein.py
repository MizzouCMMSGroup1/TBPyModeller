"""

Proein class for TBModeller

Copyright (c) Sean Lander 2014 under GPL v2

"""

from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import Entrez

import sqlite3

class Protein:

	# Structure -> Model -> Chain -> Residue -> Atom

	_pid = None
	_seq = None
	_pdb = None
	_templates = []

	def __init__(self,pid=None,sequence=None):
		self.pid = pid
		self.seq = sequence
		self.pdb = {"REMARK":[],"ATOMS":[]}
		if self.pid:
			self.pdb["REMARK"].append("Template for target {0}".format(self.pid))

	def saveToDb(self):
		print("Save to db")

	def getTemplates(self):
		print(self.seq)
		result_handle = NCBIWWW.qblast("blastp","pdb",str(self.seq),expect=0.01)
		# http://biopython.org/DIST/docs/api/Bio.Blast.NCBIWWW-module.html
		blast_records = NCBIXML.parse(result_handle)
		for record in blast_records:
			for alignment in record.alignments:
				print("SequenceID:"+alignment.accession)
				for hsp in alignment.hsps:
					print('****ALIGNMENT***')
					print('sequence:', alignment.title)
					print('length:',alignment.length)
					print('e value:', hsp.expect)
					print('score:', hsp.score)
					print('identities:',(100*hsp.identities/len(self.seq))) # need to print percentage of similarities
					print("Target:" + hsp.query[0:75] + '...')
					print("Match:" + hsp.match[0:75] + '...')
					print(hsp.sbjct[0:75] + '...')

	def choosePDB(self):
		self.pdb = self.templates[0]

	def __str__(self):
		return str(self.pdb)
