from modeller import *
from modeller.automodel import *

env = environ()
env.io.atom_files_directory = ['pdb_files']


a = automodel(env, alnfile='T0644-3U1WA.ali',
              knowns='3U1WA', sequence='T0644',
              assess_methods=(assess.DOPE, assess.GA341))
a.starting_model = 1
a.ending_model = 5
a.make()