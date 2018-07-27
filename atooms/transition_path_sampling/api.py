import numpy as np
from atooms.trajectory.decorators import Unfolded
from atooms.simulation import Simulation
from atooms.system import Thermostat
from atooms.backends.lammps import LAMMPS
from atooms.core.utils import setup_logging
from atooms.transition_path_sampling import core, TransitionPathSampling
from atooms.trajectory import TrajectoryXYZ


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
    return K

def write_thermo_tps(sim):
    """Write basic thermodynamic data."""
    f = sim.output_path + '.thermo'
    if sim.current_step == 0:
        with open(f, 'w') as fh:
            fh.write('# columns: steps, order parameter, normalized order parameter\n')
    else:
        q = sim.sim[0].order_parameter
        #q = core.calculate_order_parameter(sim.trj[i])
        with open(f, 'a') as fh:
            fh.write('%d %g %g\n' % (sim.current_step, q, q / len(sim.sim[0].system.particle)))

def main(output, input_file=None, field=0.0, steps=0, T=-1.0,
         dt=0.005, frames=60, delta_t=-1.0, t_obs=-1.0, script='',
         verbose=False):

    if input_file is None:
        raise ValueError('Provide input file')

    if T <= 0:
        raise ValueError('Provide temperature')

    if verbose:
        setup_logging(name='atooms.simulation', level=20)
        setup_logging(name='transition_path_sampling', level=40)

    nsim = 1
    if delta_t > 0:
        t_obs = delta_t * frames
    elif t_obs > 0:
        delta_t = t_obs / frames
    else:
        raise ValueError('Provide delta_t or t_obs')

    # Local dictionary to interpolate output path
    _db = locals()
    for key in _db:
        print('# {}: {}'.format(key, _db[key]))

    # Source lammps command
    import os
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
            #fix             1 all nvt temp {0} {0} 100.0
            fix             1 all nve
            timestep        {1}
            """.format(T, dt)
        else:
            cmd = script
        
    # Prepare backends
    sim= []
    for i in range(nsim):
        lmp = LAMMPS(input_file, cmd)
        lmp.verbose = False
        sim.append(Simulation(lmp, steps=int(round(delta_t / dt))))

    # Interpolate output path with input parameters
    # Ex.: output = 'output_s{field}_tobs{t_obs}'
    output = output.format(**_db)

    # Setup and run TPS simulation
    tps = TransitionPathSampling(sim, output_path=output, temperature=T, steps=steps, frames=frames, biasing_field=field)
    for s in tps.sim:
        s.system.thermostat = Thermostat(T)
    tps.add(write_thermo_tps, 1)
    core.calculate_order_parameter = mobility
    tps.run()

