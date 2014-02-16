from modeller.automodel import *
log.verbose()
env = environ()

#env.io.atom_files_directory = 'pdb_files'
env.io.atom_files_directory = 'pdb2'

a = automodel( env, alnfile='build_profile.ali', knowns='3U1W', sequence='T0644')

a.starting_model=1
a.ending_model=1

a.make()