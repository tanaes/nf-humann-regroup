#!/usr/bin/env python

import sys
import scipy as sp
from biom import Table
from biom.util import biom_open

args = sys.argv
mtx_fp = args[1]
obs_fp = args[2]
ids_fp = args[3]
out_fp = args[4]

with open(ids_fp, 'r') as f:
    ids = [x.rstrip() for x in f.readlines()]

with open(obs_fp, 'r') as f:
    obs = [x.rstrip() for x in f.readlines()]

mtx = sp.io.mmread(mtx_fp)

tab = Table(mtx, rows, cols)

with biom_open(out_fp, 'w') as f:
    tab.to_hdf5(f, 'CuratedMetagenomicData')

