# Import modules
import os
import sys
import numpy
import re
import shutil
import glob
import math
from aux_seqbuild import def_vals # function definitions
#------------------------------------------------------------------

# Set defaults
graft_opt = []; backbone_seq = []
input_pdb = 'none'; input_namd = 'none'; input_prm = 'none'
itertype = 'single'
def_res = 'none'; seg_name = 'SEG'; res_initiator = 'none'
casenum,mono_deg_poly,num_chains,fpdbflag,ftopflag,disperflag,\
fresflag,makepdifile,fnamdflag,pmolflag,cleanslate,packtol = def_vals()

print(casenum, mono_deg_poly,cleanslate,packtol)
