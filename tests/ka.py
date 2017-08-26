#!/usr/bin/env python

"""Test TPS simulation for KA mixture."""

import numpy as np
from copy import copy
from atooms.trajectory import TrajectoryRam, TrajectoryXYZ
from atooms.trajectory.decorators import Unfolded
from atooms.simulation import Simulation
from atooms.simulation.backend_lammps import LammpsBackend
from atooms.utils import setup_logging
import tps

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

# TPS parameters
tps.runsteps = 10000
tps.temperature = 0.8
tps.k = 0.01
tps.calculate_order_parameter = mobility

# Parameters
n = 1
slices = 2
steps = 10000
file_inp = 'data/ka_rho1.2.xyz'
setup_logging(level=10)  # 20 is verbose, 40 just warnings

# Set up data structures
# Replicas of simulations (each simulation might have different parameters)
cmd = """
pair_style      lj/cut 2.5
pair_coeff      1 1 1.0 1.0  2.5
pair_coeff      1 2 1.5 0.8  2.0
pair_coeff      2 2 0.5 0.88 2.2
neighbor        0.3 bin
neigh_modify    every 20 delay 0 check no
fix             1 all nve
"""
sim = [Simulation(LammpsBackend(file_inp, cmd)) for i in range(n)]
# Trajectories objects, one per simulation instance. They will have 0 frames each.
trj = [TrajectoryRam() for i in range(len(sim))]
# Input trajectories, we only read the first frame to initialize the systems
inp = [TrajectoryXYZ(file_inp) for i in range(len(sim))]

# values of the order parameter for every biasing umbrella
umbrellas = range(len(sim))

# initial value of the pseudopotential for every initial trajectory
# ideally, it should be a variable of the trajectory class
bias = np.zeros(len(sim))
# FT initialise trajectories: is there a more elegant way?
for i in range(len(sim)):
    #frame 0:
    sim[i].system = copy(inp[i][0])  # copy() might not be necessary here
    sim[i].run(steps)
    # Trajectory frame assignement takes care of copying, so copy() is not necessary here
    trj[i][0] = sim[i].system
    # Run 0
    for j in range(1, slices):
        sim[i].system = copy(trj[i][j-1])  # copy() might not be necessary here
        sim[i].run(steps)
        trj[i][j] = sim[i].system

# TPS simulation
tps_steps = 20
# calculating the value of the initial potential
tps.setup(bias, sim,trj,umbrellas)

print "bias after initialisation",bias
for k in range(tps_steps):
    tps._log.info('tps step %s', k)
    # We might have several replicas of simulations with different parameters
    for i in range(len(sim)): #FT: to be distributed?
        tps.mc_step(sim[i], trj[i], umbrellas[i], bias[i])


