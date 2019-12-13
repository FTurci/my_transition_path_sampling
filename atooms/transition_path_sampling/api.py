import os
import logging
import numpy as np
import atooms.core.progress
from atooms.trajectory.decorators import Unfolded, filter_species
from atooms.simulation import Simulation
from atooms.system import Thermostat
from atooms.backends.lammps import LAMMPS
from atooms.core.utils import setup_logging, mkdir
from atooms.transition_path_sampling import core, TransitionPathSampling
from atooms.trajectory import TrajectoryXYZ


log = logging.getLogger(__name__)


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
    unfoldedtj.add_callback(filter_species, '1')
    pos_0 = unfoldedtj[0].dump('pos')
    for j in range(1, len(t)):
        pos_1 = unfoldedtj[j].dump('pos')
        K += np.sum((pos_1 - pos_0)**2)
        pos_0 = pos_1
    t.callbacks.pop()
    return K

def write_thermo_tps(sim):
    """Write basic thermodynamic data."""
    f = sim.output_path + '.thermo'
    if sim.current_step == 0:
        with open(f, 'w') as fh:
            fh.write('# columns: steps, order parameter q, q / N, q / (N*t_obs)\n')
    else:
        q = sim.sim[0].order_parameter
        with open(f, 'a') as fh:
            fh.write('%d %g %g %g\n' % (sim.current_step, q, q / (len(sim.sim[0].system.particle)),
                                        q / (len(sim.sim[0].system.particle) * sim._tobs)))

def write_msd_tps(sim):
    f = sim.output_path + '.xyz'
    if sim.current_step == 0:
        with open(f, 'w') as fh:
            pass
    with open(f, 'a') as fh:
        fh.write('# step = %d\n' % (sim.current_step))
        for s in Unfolded(sim.trj[0], fixed_cm=True):
            fh.write('%g %g\n' % (s.particle[0].position[0], s.particle[0].velocity[0]))
        fh.write('\n\n')

def main(output, input_file=None, field=0.0, steps=0, T=-1.0,
         dt=0.005, frames=-1, delta_t=-1.0, t_obs=-1.0, script='',
         verbose=False, shift_weight=1, shoot_weight=1, debug=False,
         trajectory_interval=0, thermo_interval=1):

    # Initial checks
    if input_file is None:
        raise ValueError('Provide input file')

    if T <= 0:
        raise ValueError('Provide temperature')

    if verbose:
        setup_logging(name='atooms.simulation', level=40)
        setup_logging(name='atooms.transition_path_sampling', level=20)
        atooms.core.progress.active = False

    if debug:
        setup_logging(name='atooms.simulation', level=40)
        setup_logging(name='atooms.transition_path_sampling', level=10)
        atooms.core.progress.active = False

    # Define time intervals
    nsim = 1
    if delta_t > 0 and frames > 0:
        t_obs = delta_t * frames
    elif t_obs > 0  and frames > 0:
        delta_t = t_obs / frames
    elif t_obs > 0  and delta_t > 0:
        frames = int(round(t_obs / delta_t))
    else:
        raise ValueError('Provide two parameters out of delta_t, t_obs, frames')

    # Local dictionary to interpolate output path
    _db = locals()

    # Interpolate output path with input parameters
    # Ex.: output = 'output_s{field}_tobs{t_obs}'
    output = output.format(**_db)
    _db['output'] = output
    mkdir(os.path.dirname(output))

    # Source lammps command
    if os.path.exists(script):
        cmd = open(script).read()
    else:
        if len(script) == 0:
            # TODO: remove this KA default
            cmd = """
            pair_style      lj/cut 2.5
            pair_coeff      1 1 1.0 1.0  2.5
            pair_coeff      1 2 1.5 0.8  2.0
            pair_coeff      2 2 0.5 0.88 2.2
            neighbor        0.3 bin
            neigh_modify    every 20 delay 0 check no
            #velocity        all create {0} 12345
            timestep        {1}
            """.format(T, dt)
        else:
            cmd = script

    # Prepare backends
    lmp = LAMMPS(input_file, cmd)
    lmp.verbose = False
    sim = Simulation(lmp, steps=int(round(delta_t / dt)))

    # Always log to file
    setup_logging(filename=output + '.log', name='atooms.transition_path_sampling', level=20)
    atooms.core.progress.active = False
    # Report local parameters db
    for key in _db:
        log.info('{:12s}: {}'.format(key, _db[key]))

    # Setup and run TPS simulation
    tps = TransitionPathSampling(sim, output_path=output,
                                 temperature=T, steps=steps,
                                 frames=frames, biasing_field=field,
                                 shift_weight=shift_weight,
                                 shoot_weight=shoot_weight)
    tps._tobs = t_obs
    tps.sim.system.thermostat = Thermostat(T, relaxation_time=10.0)
    if thermo_interval > 0:
        tps.add(write_thermo_tps, thermo_interval)
    if trajectory_interval > 0:
        tps.add(write_msd_tps, trajectory_interval)
    core.calculate_order_parameter = mobility
    tps.run()
