"""

Proein class for TBModeller

Copyright (c) Sean Lander 2014 under GPL v3

"""

import sqlite3

class Protein:

	__template = []
	__seq = None

	def __init__(self,sequence=None):
		self.seq = sequence

	def saveToDb(self):
		print("Save to db")

	def __str__(self):
		return str(self.seq)
