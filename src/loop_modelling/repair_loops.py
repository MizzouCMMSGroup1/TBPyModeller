#currently, this needs to be built by hand
#using the residue file and linked to the target
#but could be automated eventually

#add 1 to each loop for now (10-20 becomes 9-21)
#to center things better

# Loop refinement of an existing model
from modeller import *
from modeller.automodel import *

log.verbose()
env = environ()

# directories for input atom files
env.io.atom_files_directory = './::../targets'

# Create a new class based on 'loopmodel' so that we can redefine
# select_loop_atoms (necessary)
class MyLoop(loopmodel):
    # This routine picks the residues to be refined by loop modeling
    def select_loop_atoms(self):
        # 10 residue insertion 
#        return selection(self.residue_range('1:', '23:'))
#        return selection(self.residue_range('121', '125'))
        return selection(self.residue_range('1:', '34:'),
                 self.residue_range('121:', '125:'))

m = MyLoop(env,
           inimodel='T0644_chain.pdb', # initial model of the target
           sequence='T0644')          # code of the target

m.loop.starting_model= 1           # index of the first loop model 
m.loop.ending_model  = 1          # index of the last loop model
m.loop.md_level = refine.slow # loop refinement method; this yields
                                   # models quickly but of low quality;
                                   # use refine.slow for better models

m.make()