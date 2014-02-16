"""

Proein class for TBModeller

Copyright (c) Sean Lander 2014 under GPL v2

"""

from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO
from Bio import Entrez
Entrez.email = "asdf@example.com"
from Bio import File

import sqlite3
import urllib
import io
import os
import time
import re
import sys
import json

rcsb_url = 'http://www.rcsb.org/pdb/download/downloadFile.do?fileFormat=pdb&compression=NO&structureId=%s'  ## PDB website

class Protein:

	# Structure -> Model -> Chain -> Residue -> Atom

	NMR_RESOLUTION = 3.5

	_debug = False
	_pid = None
	_seq = None
	_pdb = None
	_pdbs = None
	_templates = None
	_templatesfolder = None
	_fastas = None
	_fastasfolder = None
	_alignments = None
	_alignmentsfolder = None

	def __init__(self,pid=None,sequence=None,debug=False):
		self.pid = pid if pid else str(datetime.now())
		self.seq = sequence
		self.debug = debug
		self.templates = {}
		self.fastas = {}
		self.alignments = {}
		self.pdb = None
		self.pdbs = []
		# Prepare folders
		self.fastasfolder = 'fastas-%s' % str(self.pid)
		self.alignmentsfolder = 'alignments-%s' % str(self.pid)
		self.templatesfolder = 'templates-%s' % str(self.pid)
		if not os.path.exists(self.fastasfolder):
			os.mkdir(self.fastasfolder)
		if not os.path.exists(self.alignmentsfolder):
			os.mkdir(self.alignmentsfolder)
		if not os.path.exists(self.templatesfolder):
			os.mkdir(self.templatesfolder)
		self.ex_resolution = re.compile(b'REMARK\s*2 RESOLUTION\.\s*([0-9\.]+|NOT APPLICABLE)\s+.*')

	def saveToDb(self):
		print("Save to db")

	def getTemplates(self):
		# http://biopython.org/DIST/docs/api/Bio.Blast.NCBIWWW-module.html
		if self.debug:
			print(self.seq)
			result_handle = NCBIWWW.qblast("blastp","pdb",str(self.seq))
			print("BLAST Request Finished")
			sys.stdout.flush()
		else:
			result_handle = NCBIWWW.qblast("blastp","pdb",str(self.seq),expect=0.01)
		blast_records = NCBIXML.parse(result_handle)
		for record in blast_records:
			for alignment in record.alignments:
				id = alignment.accession
				self.templates[id] = {"fasta":self.getFastaFromId(id),
					'asequence':alignment.title,
					'alength':alignment.length,
					"alignments":[]}
				for hsp in alignment.hsps:
					self.templates[id]["alignments"].append({'expect':hsp.expect,
						'score':hsp.score,
						'identities':hsp.identities,
						'similarity':(100*hsp.identities/len(self.seq)),
						'target':hsp.query,
						'match':hsp.match,
						'template':hsp.sbjct})
					if self.debug:
						print()
						print('****ALIGNMENT***')
						print('id:',alignment.accession)
						print('sequence:', alignment.title)
						print('length:',alignment.length)
						print('e value:', hsp.expect)
						print('score:', hsp.score)
						print('identities:',(100*hsp.identities/len(self.seq))) # need to print percentage of similarities
						print("Target  :" + hsp.query[0:75] + '...')
						print("Match   :" + hsp.match[0:75] + '...')
						print("Template:" + hsp.sbjct[0:75] + '...')
						print()

				self.fastas[id] = self.templates[id]['fasta']
				self.alignments[id] = self.templates[id]['alignments']
				fname = '%s/%s.fasta' % (self.fastasfolder,id)
				if not os.path.exists(fname):
					f = open(fname,'w')
					SeqIO.write(self.fastas[id],f,'fasta')
					f.close()
				for i,a in enumerate(self.alignments[id]):
					fname = '%s/%s-%s.alignment' % (self.alignmentsfolder,id,str(i))
					if not os.path.exists(fname):
						f = open(fname,'w')
						json.dump(a,f)
						f.close()

	def getFastaFromId(self,id):
		handle = Entrez.efetch(db="protein", rettype="fasta", id=id)
		frecord = SeqIO.read(handle, "fasta")
		frecord.id = str(id)
		if self.debug:
			print(frecord)
		handle.close()
		return frecord

	def getPDBs(self):
		results = []
		for id in self.templates.keys():
			handle = self.getRemotePDBHandle(id.split('_')[0])
			lines,infos = self.parsePdbFromHandle(handle)
			if id in self.templates:
				self.templates[id]['fasta'].__dict__.update(infos)
			fname = '%s/%s.pdb' % (self.templatesfolder,id)
			if not os.path.exists(fname):
				f = open(fname,'wb',1)
				f.writelines(lines)
				f.close()
			results.append(fname)

	def getRemotePDBHandle( self, id ):
		"""
		From Biskit

		Get the coordinate file remotely from the RCSB.

		@param id: pdb code, 4 characters
		@type  id: str
		@param rcsb_url: template url for pdb download
		(default: L{settings.rcsb_url})
		@type  rcsb_url: str

		@return: the requested pdb file as a file handle
		@rtype: open file handle

		@raise BlastError: if couldn't retrieve PDB file
		"""
		handle = urllib.request.urlopen( rcsb_url% (id) )

		uhandle = File.UndoHandle(handle)

		if not uhandle.peekline():
			raise BaseException( "Couldn't retrieve ", rcsb_url )

		return uhandle

	def parsePdbFromHandle(self, handle, first_model_only=True ):
		"""
		From Biskit

		Parse PDB from file/socket or string handle into memory.

		@param handle: fresh open file/socket handle to PDB ressource or string
		@type  handle: open file-like object or str
		@param first_model_only: only take first of many NMR models [True]
		@type  first_model_only: bool

		@return: pdb file as list of strings, dictionary with resolution
		@rtype: [str], {'resolution':float }
		@raise BlastError: if passed in string is too short
		"""
		lines = []
		res_match = None
		infos = {}

		if type( handle ) is str:
			if len(handle) < 5000:
				raise BaseException( "Couldn't extract PDB Info." )
			handle =  io.StringIO( handle )

		## if handle.peekline()[:6] != 'TITLE':
		##     raise BlastError, 'Ressource does not seem to be a PDB:\n%r' %\
		##   handle.peekline()

		for l in handle:
			lines += [ l ]

			res_match = res_match or self.ex_resolution.search( l )

			if first_model_only and l[:6] == 'ENDMDL':
				break

		if res_match:
			if res_match.groups()[0] == 'NOT APPLICABLE':
				infos['resolution'] = self.NMR_RESOLUTION
			else:
				infos['resolution'] = float( res_match.groups()[0] )
		else:
			#raise BaseException('No resolution record found in PDB.')
			print('No resolution record found in PDB.')

		return lines, infos

	def choosePDB(self):
		self.pdb = None

	def __str__(self):
		return str(self.pid)
