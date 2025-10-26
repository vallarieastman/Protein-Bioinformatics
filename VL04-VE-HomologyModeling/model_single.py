# -*- coding: utf-8 -*-
from modeller import *
from modeller.automodel import *

log.verbose()
env = Environ()
env.io.atom_files_directory = ['.']

# Define the automodel class
a = AutoModel(env,
              alnfile='TvLDH-4UUM.ali',      # alignment from step 1
              knowns='4UUM_A',               # template code in the alignment
              sequence='TvLDH')              # target code in the alignment

a.starting_model = 1                         # how many models to build
a.ending_model  = 10                         # build 10 models
a.assess_methods = (assess.DOPE, assess.GA341)

a.make()