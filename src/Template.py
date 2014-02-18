"""

Template class for TBModeller

Copyright (c) Sean Lander 2014 under GPL v2

"""

from ModObj import ModObj
from Alignment import Alignment
from PDB import PDB

class Template(ModObj):

	'''
	id
	fasta
	sequence
	length
	alignments
	pdb
	'''

	def __init__(self,id=None,fasta=None,sequence=None,length=0,alignments=[],pdb=None):
		self.id = id
		self.fasta = fasta
		self.sequence = sequence
		self.length = length
		self.alignments = alignments
		self.pdb = pdb