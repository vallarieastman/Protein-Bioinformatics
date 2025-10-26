# -*- coding: utf-8 -*-

from modeller import *
from modeller.automodel import *

log.verbose()
env = Environ()
env.io.atom_files_directory = ['.']

aln = Alignment(env)
aln.append(file='TvLDH-4UUM.ali', alignment_format='PIR')  # optional if you generated it
# OR regenerate from scratch:
# aln = Alignment(env)
# aln.append(file='TvLDH.pir', alignment_format='PIR', align_codes='TvLDH')
# mdl = Model(env, file='4UUM', model_segment=('FIRST:A','LAST:A'))
# aln.append_model(mdl, align_codes='4UUM_A', atom_files='4UUM.pdb')
# aln.align2d()
# aln.write(file='TvLDH-4UUM.ali', alignment_format='PIR')

# proceed to modeling using the .ali you now trust
