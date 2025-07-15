# singularity run -B /nfs,/lustre /nfs/cellgeni/singularity/images/snapatac2_v2.8.0.sif python

import argparse
import snapatac2 as snap
import pandas as pd
import scanpy as sc
import numpy as np

n_jobs = 10

data = snap.read_dataset('work/results_snapatac2_combine/full_adatas/full.h5ads')
snap.metrics.tsse(data,snap.genome.hg38,n_jobs=n_jobs)
#snap.metrics.frip(data,n_jobs=n_jobs)
data.close()

