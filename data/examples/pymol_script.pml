# script: 3C_loopmap.pml
#
#setup global parameters

@pymol_setup.mac

load 4FR9.pdb, 3C

set cartoon_fancy_helices,1

cartoon loop, 3C
set cartoon_loop_radius,0.3, 3C colour tv_yellow, 3C

show cartoon, 3C

