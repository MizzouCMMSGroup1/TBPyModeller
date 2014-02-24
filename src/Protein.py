"""

Protein class for TBModeller

Copyright (c) Sean Lander 2014 under GPL v2

"""

from Template import Template
from Alignment import Alignment
from PDB import PDB
from Atom import Atom

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
import copy

rcsb_url = 'http://www.rcsb.org/pdb/download/downloadFile.do?fileFormat=pdb&compression=NO&structureId=%s'  ## PDB website

amino3to1 = {
	"ALA":"A",
	"ARG":"R",
	"ASN":"N",
	"ASP":"D",
	"CYS":"C",
	"GLU":"E",
	"GLN":"Q",
	"GLY":"G",
	"HIS":"H",
	"ILE":"I",
	"LEU":"L",
	"LYS":"K",
	"MET":"M",
	"PHE":"F",
	"PRO":"P",
	"SER":"S",
	"THR":"T",
	"TRP":"W",
	"TYR":"Y",
	"VAL":"V",
	"MYL":"X",
	"---":"-"
}
amino1to3 = {
	"A":"ALA",
	"R":"ARG",
	"N":"ASN",
	"D":"ASP",
	"C":"CYS",
	"E":"GLU",
	"Q":"GLN",
	"G":"GLY",
	"H":"HIS",
	"I":"ILE",
	"L":"LEU",
	"K":"LYS",
	"M":"MET",
	"F":"PHE",
	"P":"PRO",
	"S":"SER",
	"T":"THR",
	"W":"TRP",
	"Y":"TYR",
	"V":"VAL",
	"X":"MYL",
	"-":"---"
}

class Protein:

	# Structure -> Model -> Chain -> Residue -> Atom

	NMR_RESOLUTION = 3.5

	debug = False
	'''
	Protein variables
	pid
	seq
	templates
	templatesfolder
	fastas
	fastasfolder
	alignments
	alignmentsfolder
	pdb
	pdbs
	targetsfolder
	'''

	def __init__(self,pid=None,sequence=None,debug=False):
		self.pid = pid if pid else str(datetime.now())
		self.seq = sequence
		self.debug = debug
		self.templates = {}
		self.fastas = {}
		self.alignments = []
		self.pdb = None
		self.pdbs = {}
		# Prepare folders
		self.fastasfolder = 'fastas-%s' % str(self.pid)
		self.alignmentsfolder = 'alignments-%s' % str(self.pid)
		self.templatesfolder = 'templates-%s' % str(self.pid)
		self.targetsfolder = 'targets'
		if not os.path.exists(self.fastasfolder):
			os.mkdir(self.fastasfolder)
		if not os.path.exists(self.alignmentsfolder):
			os.mkdir(self.alignmentsfolder)
		if not os.path.exists(self.templatesfolder):
			os.mkdir(self.templatesfolder)
		if not os.path.exists(self.targetsfolder):
			os.mkdir(self.targetsfolder)
		# Resolution regex for PDB file
		self.ex_resolution = re.compile(b'REMARK\s*2 RESOLUTION\.\s*([0-9\.]+|NOT APPLICABLE)\s+.*')

	def saveToDb(self):
		print("Save to db")

	def getTemplates(self):
		# http://biopython.org/DIST/docs/api/Bio.Blast.NCBIWWW-module.html
		if self.debug:
			print(self.seq)
		# Send BLAST request to server
		# Use blastp (protein) for the method
		# Use pdb as the database
		result_handle = NCBIWWW.qblast("blastp","pdb",str(self.seq),expect=0.01)
		# Parse the results into blast records
		blast_records = NCBIXML.parse(result_handle)
		if self.debug:
			print("BLAST Request Finished")

		# Read through each blast record
		for record in blast_records:
			# Grab the alignments from each record
			for alignment in record.alignments:
				# Use the alignment id as the template key
				id = alignment.accession
				fasta = self.getFastaFromId(id)
				title = alignment.title
				length = alignment.length
				# Set up the template object for this id
				template = Template(
					id=id,fasta=fasta,sequence=title,
					length=length,alignments=[]
				)
				# Store the template in the template dict
				self.templates[id] = template
				"""
				self.templates[id] = {"fasta":self.getFastaFromId(id),
					'asequence':alignment.title,
					'alength':alignment.length,
					"alignments":[]}
				"""
				# Store fasta in dict
				self.fastas[id] = fasta
				# Get all alignments for this template
				for hsp in alignment.hsps:
					# Create an alignment object
					a = Alignment(
						id=id,title=title,expect=hsp.expect,score=hsp.score,
						identities=hsp.identities,similarity=(100*hsp.identities/len(self.seq)),
						target=hsp.query,targetstart=hsp.query_start,match=hsp.match,
						template=hsp.sbjct,templatestart=hsp.sbjct_start,length=length
					)
					# Alignment isn't necessarily the same size as the sequence
					targetfront = str(self.seq[:a.targetstart-1])
					targetend = str(self.seq[(a.targetstart+a.length):])
					a.target = ''.join(targetfront) + a.target + ''.join(targetend)
					a.length = len(a.target)
					
					templatefront = ['-']*(a.targetstart-1)
					templateend = ['-']*(len(self.seq)-(a.targetstart+a.length))
					a.template = ''.join(templatefront) + a.template + ''.join(templateend)

					if self.debug:
						print("Seq vs Target Length:",len(self.seq),a.length)

					# Append the alignment to the template's alignments
					self.templates[id].alignments.append(a)
					self.alignments.append(a)
					"""
					self.templates[id]["alignments"].append({'expect':hsp.expect,
						'score':hsp.score,
						'identities':hsp.identities,
						'similarity':(100*hsp.identities/len(self.seq)),
						'target':hsp.query,
						'match':hsp.match,
						'template':hsp.sbjct})
					"""

					if self.debug:
						print()
						print('****ALIGNMENT***')
						print('id:',id)
						print('sequence:', title)
						print('length:',length)
						print('e value:', hsp.expect)
						print('score:', hsp.score)
						print('identities:',(100*hsp.identities/len(self.seq))) # need to print percentage of similarities
						print("Target  :" + hsp.query[0:75] + '...')
						print("Match   :" + hsp.match[0:75] + '...')
						print("Template:" + hsp.sbjct[0:75] + '...')
						print()

		# Save off the fasta file
		for id,fasta in self.fastas.items():
			fname = '%s/%s.fasta' % (self.fastasfolder,id)
			if not os.path.exists(fname):
				f = open(fname,'w')
				SeqIO.write(fasta,f,'fasta')
				f.close()

		# Save off the alignments
		for i,a in enumerate(self.alignments):
			fname = '%s/%s-%s.alignment' % (self.alignmentsfolder,a.id,str(i))
			if not os.path.exists(fname):
				f = open(fname,'w')
				json.dump(a.toJSON(),f)
				f.close()

		return self.templates.keys()

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
			pdb = PDB(infos=infos,lines=lines)
			if id in self.templates:
				self.templates[id].fasta.__dict__.update(infos)
			self.templates[id].pdb = pdb
			self.pdbs[id] = pdb
			fname = '%s/%s.pdb' % (self.templatesfolder,id)
			if not os.path.exists(fname):
				f = open(fname,'wb',1)
				f.writelines(lines)
				f.close()
			results.append(fname)
		return results

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

			if first_model_only and l[:6] == b'ENDMDL':
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

	def alignPDB(self, skipBlanks = True):
		# Sort the alignments by similarity
		alignments = self.alignments
		alignments.sort(key = lambda a: a.score,reverse=True)
		target = None
		template = None
		alignment = None
		# Get the highest scoring alignment with similarity below 90%
		for a in alignments:
			if a.similarity > 90: continue
			else:
				alignment = a
				target = a.target
				template = a.template
				pdb = self.pdbs[a.id]
				break
		if self.debug:
			print("Selected Template",alignment.id)
		# We now have the alignment and the pdb
		# Convert the pdb to a residue sequence
		# DO NOT CHANGE THESE
		seq = []
		curres = '---'
		res = []
		first = True
		atomnum = 5
		resnum = 0
		for line in pdb.lines:
			# Decode from byte array
			line = line.decode()
			# Grab only ATOM lines
			if line[:6] == 'ATOM  ':
				# Get the atomname (N,CA,C,O)
				atomname = line[12:16]

				# First make sure you're only storing N, CA, C and O
				if atomname not in (' N  ',' CA ',' C  ',' O  ') or line[16] not in (' ','A'):
					continue

				# Check to see if residue changed
				if line[17:20] != curres or atomnum > 4:
					if self.debug:
						print("%8s%8s%8s%8s" % ('CURRES','NEWRES','RESNUM','ATONUM'))
						print("%8s%8s%8d%8d" % (curres,line[17:20],resnum,atomnum))
						print(line)
					# Check to make sure the last residue was complete
					if atomnum < 5:
						# Add on the missing atoms to the end of the residue
						if atomnum == 1:
							atom = Atom(atomName=' N  ',resName=curres,elemSym=' N',missing=True)
							res.append(atom)
							atomnum += 1
						if atomnum == 2:
							atom = Atom(atomName=' CA ',resName=curres,elemSym=' C',missing=True)
							res.append(atom)
							atomnum += 1
						if atomnum == 3:
							atom = Atom(atomName=' C  ',resName=curres,elemSym=' C',missing=True)
							res.append(atom)
							atomnum += 1
						if atomnum == 4:
							atom = Atom(atomName=' O  ',resName=curres,elemSym=' O',missing=True)
							res.append(atom)
							atomnum += 1
					# Append the residue to the sequence
					if res: seq.append((amino3to1[curres],res))

					# If this is the first atom fill in the blanks from the PDB
					if first:
						first = False
						if self.debug:
							print("Starting resnum:",int(line[22:26]))
						# Residue number of first residue in PDB
						for k in range(1,int(line[22:26])):
							resname = '---'
							res = []
							res.append(Atom(atomName=' N  ',resName=resname,elemSym=' N',missing=True))
							res.append(Atom(atomName=' CA ',resName=resname,elemSym=' C',missing=True))
							res.append(Atom(atomName=' C  ',resName=resname,elemSym=' C',missing=True))
							res.append(Atom(atomName=' O  ',resName=resname,elemSym=' O',missing=True))
							seq.append((amino3to1[resname],res))
						# Decrement resnum by one to mimic the last res not in the PDB
						resnum = int(line[22:26])-1

					# Check for resnum increment so we don't miss residues
					resinc = int(line[22:26])-resnum
					if resinc != 1:
						if self.debug:
							print("Resiude skip:",resnum,"->",int(line[22:26]))
						for k in range(1,resinc):
							resname = '---'
							res = []
							res.append(Atom(atomName=' N  ',resName=resname,elemSym=' N',missing=True))
							res.append(Atom(atomName=' CA ',resName=resname,elemSym=' C',missing=True))
							res.append(Atom(atomName=' C  ',resName=resname,elemSym=' C',missing=True))
							res.append(Atom(atomName=' O  ',resName=resname,elemSym=' O',missing=True))
							seq.append((amino3to1[resname],res))

					# Reset the residue
					res = []
					# Store residue numbers
					resnum = int(line[22:26])
					# Store the atom number
					atomnum = 1

				# Make sure current residue name is set
				curres = line[17:20]

				'''
				if self.debug:
					print("Line:",line)
				'''

				# Check to make sure the right atom is being stored
				if atomnum == 1:
					if atomname != ' N  ':
						atom = Atom(atomName=' N  ',resName=curres,elemSym=' N',missing=True)
						res.append(atom)
						atomnum += 1
				if atomnum == 2:
					if atomname != ' CA ':
						atom = Atom(atomName=' CA ',resName=curres,elemSym=' C',missing=True)
						res.append(atom)
						atomnum += 1
				if atomnum == 3:
					if atomname != ' C  ':
						atom = Atom(atomName=' C  ',resName=curres,elemSym=' C',missing=True)
						res.append(atom)
						atomnum += 1
				if atomnum == 4:
					if atomname != ' O  ':
						atom = Atom(atomName=' O  ',resName=curres,elemSym=' O',missing=True)
						res.append(atom)
						atomnum += 1
						continue
				if atomnum > 4:
					continue
				# Store off the atom information if valid
				atom = Atom(
					atomName=atomname,
					altLoc=line[16],
					resName=line[17:20],
					chainId=line[21],
					codeForInsertion=line[26],
					xcoord=float(line[30:38]),
					ycoord=float(line[38:46]),
					zcoord=float(line[46:54]),
					occ=float(line[54:60]),
					temp=float(line[60:66]),
					elemSym=line[76:78],
					charge=line[78:80]
				)
				# And add it to the residue (in the proper order)
				res.append(atom)
				atomnum += 1

		# Add the final residue to the sequence
		seq.append((amino3to1[curres],res))

		# Sequence has been pulled out
		if self.debug:
			output = []
			for s in seq:
				output.append(s[0])
			print(''.join(output))

		i = 0 # template
		j = 0 # seq
		targetseq = [] # final sequence
		if self.debug:
			print("%4s%4s%10s%10s%10s" % ("i","j","TARGET","TEMPLATE","PDB"))
		while i < len(template) and j < len(seq):
			if self.debug:
				print("%4d%4d%10s%10s%10s" % (i,j,target[i],template[i],seq[j][0]))

			# Template has a blank
			# Fill in blank data for that residue
			if '-' == template[i]:
				resname = amino1to3[target[i]]
				res = []
				res.append(Atom(atomName=' N  ',resName=resname,elemSym=' N',missing=True))
				res.append(Atom(atomName=' CA ',resName=resname,elemSym=' C',missing=True))
				res.append(Atom(atomName=' C  ',resName=resname,elemSym=' C',missing=True))
				res.append(Atom(atomName=' O  ',resName=resname,elemSym=' O',missing=True))
				targetseq.append(res)
				i += 1
				if '-' == seq[j][0]:
					j += 1
				continue
			# Target has a blank
			# Skip that line
			if '-' == target[i]:
				i += 1
				j += 1
				continue
			# PDB doesn't match template sequence ERROR
			if seq[j][0] != '-' and template[i] != seq[j][0]:
				print("Template doesn't match PDB: (%04d)%c - (%04d)%c" % (i,template[i],j,seq[j][0]))

			# Everything lines up!
			# Copy over the data
			resname = amino1to3[target[i]]
			res = []
			# Take the atoms from the residue and swap residue names
			for atom in seq[j][1]:
				if not atom.missing:
					a = copy.copy(atom)
					a.resName = resname
					res.append(a)
			# Append it to the list and increment counters
			targetseq.append(res)
			i += 1
			j += 1

		# Time to make the PDB
		lines = []
		lines.append("REMARK Template for target %s\n" % (self.pid))
		residues_needed_for_loop_modelling = []
		residues_needed_for_loop_modelling.append("Loops needed for target %s\n" % (self.pid))
        
		i_atom = 0
		i_residue = 0
		last_residue_for_loop_modelling = 0
		for res in targetseq:
			i_residue += 1
			for atom in res:
				i_atom += 1
				if atom.missing:
					if (last_residue_for_loop_modelling!=i_residue):
						residues_needed_for_loop_modelling.append("%i, " % i_residue)
						last_residue_for_loop_modelling = i_residue
				# If skipBlanks is True then only print atoms that aren't missing
				# Otherwise print all atoms
				if not skipBlanks or not atom.missing:
					lines.append('ATOM  %5d %4s%c%3s %c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s\n' % (
						i_atom,atom.atomName,atom.altLoc,atom.resName,' ',i_residue,atom.codeForInsertion,
						atom.xcoord,atom.ycoord,atom.zcoord,atom.occ,atom.temp,atom.elemSym,atom.charge
					))
					 '''blank character replaces atom.chainId'''

		fname = '%s/%s.pdb' % (self.targetsfolder,self.pid)
		f = open(fname,'w')
		f.writelines(lines)
		f.close()
        
		resname = '%s/%s.residue' % (self.targetsfolder,self.pid)
		f = open(resname,'w')
		f.writelines(residues_needed_for_loop_modelling)
		f.close()

	def __str__(self):
		return str(self.pid)