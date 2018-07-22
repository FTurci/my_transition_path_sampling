#!/usr/bin/env python

"""Test TPS simulation for KA mixture."""

import numpy as np
from atooms.trajectory.decorators import Unfolded
from atooms.simulation import Simulation, target_rmsd, write_thermo, write_config
from atooms.system import Thermostat
from atooms.backends.lammps import LAMMPS
from atooms.core.utils import setup_logging
from atooms.transition_path_sampling import core, TransitionPathSampling
from atooms.trajectory import TrajectoryXYZ


setup_logging(name='atooms.simulation', level=40)  # 20 is verbose, 40 just warnings
setup_logging(name='transition_path_sampling', level=40)


def self_overlap(r0, r1, side, a_square):
    rij = np.sum((r0 - r1)**2, axis=1) # square displacement
    qij = rij.flatten() < a_square
    return qij.sum()


def mobility(t):
    # Note: the trajectory must be unfolded
    # pos = [s.dump('pos') for s in ]
    side = t[0].cell.side
    K = 0
    unfoldedtj = Unfolded(t, fixed_cm=True)
    pos_0 = unfoldedtj[0].dump('pos')
    for j in range(1, len(t)):
        pos_1 = unfoldedtj[j].dump('pos')
        K += np.sum((pos_1 - pos_0)**2)
        pos_0 = pos_1

    print ('# tps K=%s [len=%s]' % (K, len(t)))
    return K

def write_thermo_tps(sim):
    """Write basic thermodynamic data."""
    f = sim.output_path + 'thermo'
    if sim.current_step == 0:
        with open(f, 'w') as fh:
            fh.write('# columns: steps, order parameter, normalized order parameter\n')
    else:
        q = sim.sim[0].order_parameter
        #q = core.calculate_order_parameter(sim.trj[i])
        with open(f, 'a') as fh:
            fh.write('%d %g %g\n' % (sim.current_step, q, q / len(sim.sim[0].system.particle)))

nsim = 1
dt = 0.005
frames = 60
delta_t = 1.5
t_obs = delta_t*frames
print "# t_obs",t_obs
print "# granularity delta_t ", delta_t
T = 0.6

file_inp = 'data/ka_rho1.2_N150_T0.6.xyz'
cmd = """
pair_style      lj/cut 2.5
pair_coeff      1 1 1.0 1.0  2.5
pair_coeff      1 2 1.5 0.8  2.0
pair_coeff      2 2 0.5 0.88 2.2
neighbor        0.3 bin
neigh_modify    every 20 delay 0 check no
#velocity        all create {0} 12345
#fix             1 all nvt temp {0} {0} 100.0
fix             1 all nve
timestep        {1}
""".format(T, dt)

import sys
field = float(sys.argv[1])
# Initial equilibration
# equiibration_rmsd = 0.1
# check_frequency = 100 #steps

# equilibrator = Simulation(LAMMPS(file_inp, cmd), steps=int(t_obs/dt),output_path="eq_") 
# equilibrator.trajectory = TrajectoryXYZ
# # get the root mean square displacement to a certain value to equilibrate
# # equilibrator.add(target_rmsd, check_frequency, equiibration_rmsd)
# equilibrator.add(write_thermo,check_frequency)
# equilibrator.add(write_config,check_frequency)
# equilibrator.run(int(10000))
# write final configuration
# with TrajectoryXYZ(file_inp+"_equilibration.xyz", 'w') as thd:
#     thd.write(equilibrator.system,step=0)

# sim = [Simulation(LAMMPS(file_inp, cmd), steps=int(round(delta_t / dt)) ) for i in range(nsim)]
sim= []
for i in range(nsim):
    lmp = LAMMPS(file_inp, cmd)
    lmp.verbose = False
    sim.append(Simulation(lmp, steps=int(round(delta_t / dt)) ) )
tps = TransitionPathSampling(sim, output_path='output.s%g.'%field, temperature=T, steps=10000, frames=frames, biasing_field=field)
for s in tps.sim:
    s.system.thermostat = Thermostat(T)
tps.add(write_thermo_tps, 1)
core.calculate_order_parameter = mobility
tps.run()


