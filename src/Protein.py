"""

Proein class for TBModeller

Copyright (c) Sean Lander 2014 under GPL v3

"""

import sqlite3

class Protein:

	_pid = None
	_seq = None
	_pdb = None

	def __init__(self,pid=None,sequence=None):
		self.pid = pid
		self.seq = sequence
		self.pdb = {"REMARK":[],"ATOM":[]}
		if self.pid:
			self.pdb["REMARK"].append("Template for target {0}".format(self.pid))

	def saveToDb(self):
		print("Save to db")

	def __str__(self):
		return str(self.pdb)
