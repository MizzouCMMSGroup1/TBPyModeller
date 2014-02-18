"""

PDB class for TBModeller

Copyright (c) Sean Lander 2014 under GPL v2

"""

from ModObj import ModObj

class PDB(ModObj):

	'''
	infos
	lines
	'''

	def __init__(self,infos={},lines=[]):
		self.infos = infos
		self.lines = lines