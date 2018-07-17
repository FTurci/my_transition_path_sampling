#!/usr/bin/env python

"""Test TPS simulation for KA mixture."""

import numpy as np
from atooms.trajectory.decorators import Unfolded
from atooms.simulation import Simulation
from atooms.backends.lammps import LAMMPS
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
    pos = [s.dump('pos') for s in Unfolded(t, fixed_cm=True)]
    side = t[0].cell.side
    a = 0.3
    s = 0
    K = 0
    for j in range(1, len(t)):
        K += np.sum((pos[j] - pos[j-1])**2)
        #s += self_overlap(pos[j], pos[j-1], side, a**2)
    print '# tps C=%s [len=%s]' % (K, len(t))
    return K

def write_thermo(sim):
    """Write basic thermodynamic data."""
    f = sim.output_path + 'thermo'
    if sim.current_step == 0:
        with open(f, 'w') as fh:
            fh.write('# \n')
    else:
        q = sim.sim[0].order_parameter = core.calculate_order_parameter(sim.trj[i])
        with open(f, 'a') as fh:
            fh.write('%d %s\n' % (sim.current_step, q))

nsim = 1
t_obs = 5.
dt = 0.004
T = 0.6
file_inp = 'data/ka_rho1.2.xyz'
cmd = """
pair_style      lj/cut 2.5
pair_coeff      1 1 1.0 1.0  2.5
pair_coeff      1 2 1.5 0.8  2.0
pair_coeff      2 2 0.5 0.88 2.2
neighbor        0.3 bin
neigh_modify    every 20 delay 0 check no
fix             1 all nvt temp {0} {0} 100.0
# fix             1 all nve
timestep        {1}
""".format(T, dt)

sim = [Simulation(LAMMPS(file_inp, cmd), steps=int(t_obs/dt)) for i in range(nsim)]
tps = TransitionPathSampling(sim, output_path='./', temperature=0.8, steps=10, slices=10, biasing_field=0.0)
tps.add(write_thermo, 1)
core.calculate_order_parameter = mobility
tps.run()


