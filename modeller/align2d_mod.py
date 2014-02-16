from modeller import *

env = environ()
aln = alignment(env)

pdb = '3U1W'
chain = 'A'

env.io.atom_files_directory = ['pdb_files']


mdl = model(env, file=pdb, model_segment=('FIRST:A','LAST:A'))
aln.append_model(mdl, align_codes='3U1WA', atom_files='3U1W.pdb')
aln.append(file='alignment.ali', align_codes='T0644')
aln.align2d()
aln.write(file='T0644-3U1WA.ali', alignment_format='PIR')
aln.write(file='T0644-3U1WA.pap', alignment_format='PAP')