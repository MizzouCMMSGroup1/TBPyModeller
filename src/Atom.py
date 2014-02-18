"""

Atom class for TBModeller

Copyright (c) Sean Lander 2014 under GPL v2

"""

from ModObj import ModObj

class Atom(ModObj):

	def __init__(self,atomName=None,altLoc=None,resName=None,chainId=None,
		codeForInsertion=None,xcoord='0000.000',ycoord='0000.000',zcoord='0000.000',
		occ='000.00',temp='000.00',elemSym=None,charge=None):
		self.atomName = atomName
		self.altLoc = altLoc
		self.resName = resName
		self.chainId = chainId
		self.codeForInsertion = codeForInsertion
		self.xcoord = xcoord
		self.ycoord = ycoord
		self.zcoord = zcoord
		self.occ = occ
		self.temp = temp
		self.elemSym = elemSym
		self.charge = charge