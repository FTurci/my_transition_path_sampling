"""Writer callbacks specific to transition path sampling."""

import os
import numpy
from atooms.utils import rank, size, comm, barrier
from atooms.utils import rmd, rmf, mkdir

def write_thermo(tps):
    for i in range(len(tps.sim)):
        f = os.path.join(tps.output_path, '%d.out' % i)
        if tps.steps == 0:
            with open(f, 'w') as fh:
                obs = ', '.join(['steps', 'bias', 'umbrella'])
                fh.write('# columns:' + obs + '\n')
        with open(f, 'a') as fh:
            fh.write('%d %s %s\n' % (tps.steps, tps.bias[i], tps.umbrellas[i]))
