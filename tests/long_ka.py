#!/usr/bin/env python

"""Test TPS simulation for KA mixture."""

import numpy as np
from atooms.trajectory.decorators import Unfolded
from atooms.simulation import Simulation
from atooms.simulation.lammps import LAMMPS
from atooms.core.utils import setup_logging
from atooms.transition_path_sampling import core, TransitionPathSampling

setup_logging(name='atooms.simulation', level=40)  # 20 is verbose, 40 just warnings
setup_logging(name='transition_path_sampling', level=20)


def self_overlap(r0, r1, side, a_square):
    rij = np.sum((r0 - r1)**2, axis=1) # square displacement
    qij = rij.flatten() < a_square
    return qij.sum()

def mobility(t):
    # Note: the trajectory must be unfolded
    pos = [s.dump('pos') for s in Unfolded(t)]
    side = t[0].cell.side
    a = 0.3
    s = 0
    for j in range(1, len(t)):
        s += self_overlap(pos[j], pos[j-1], side, a**2)
    print '# tps C=%s [len=%s]' % (s, len(t))
    return s

nsim = 1
file_inp = 'data/ka_rho1.2.xyz'
cmd = """
pair_style      lj/cut 2.5
pair_coeff      1 1 1.0 1.0  2.5
pair_coeff      1 2 1.5 0.8  2.0
pair_coeff      2 2 0.5 0.88 2.2
neighbor        0.3 bin
neigh_modify    every 20 delay 0 check no
fix             1 all nve
"""
sim = [Simulation(LAMMPS(file_inp, cmd), steps=100) for i in range(nsim)]
tps = TransitionPathSampling(sim, temperature=0.8, steps=10, slices=10)
core.calculate_order_parameter = mobility
tps.run()


