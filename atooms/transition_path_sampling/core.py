"""Sketch of TPS simulation."""

import logging
import numpy as np
from copy import copy
from atooms.trajectory import TrajectoryRam, TrajectoryXYZ
from atooms.simulation import Simulation
from atooms.backends.dryrun import DryRun
from atooms.transition_path_sampling import __version__, __date__, __commit__

log = logging.getLogger(__name__)


def firstHalf(tj, margin=0):
    """Return an integer in [margin,trajectoryLength/2)"""
    trajectoryLength = len(tj)
    return margin+np.random.randint(trajectoryLength/2-margin)

def secondHalf(tj, margin=0):
    """Return an integer in [trajectoryLength/2,trajectoryLength-margin)"""
    trajectoryLength = len(tj)
    return trajectoryLength-margin-1-np.random.randint(trajectoryLength/2-margin)

def shootForward(sim, tj, slice):
    """
    Perform a forward shooting move: from a selected slice, randomise
    the velocities and integrate a new piece of trajectory forward in
    time.
    """
    log.debug('tps shoot forward')
    sim.system = tj[slice]
    sim.system.temperature = sim.thermostat_temperature

    # overwrtiting from slice till the end
    for j in range(slice, len(tj)):
        sim.run()
        tj[j] = sim.system
    return tj

def shootBackward(sim, tj, slice):
    """
    Perform a backward shooting move: from a selected slice, randomise
    the velocities and integrate a new piece of trajectory backward in
    time.
    """
    log.debug('tps shoot backward')
    sim.system = tj[slice]
    sim.system.temperature = sim.thermostat_temperature
    # !!!
    # WARNING: TIME REVERSAL MISSING
    #can/should i reverse time?
    for j in range(slice-1, -1, -1):
        sim.run()
        tj[j] = sim.system
    return tj

def shiftForward(sim, tj, slice):
    """
    Perform a forward shifting move: delete a piece of trajectory
    starting from slice reverse the order of the trajectory and
    continue the trajectory for the missing bit on the other end.
    """
    log.debug('tps shift forward')
    trajectoryLength = len(tj)
    copytj = TrajectoryRam()
    # reverse copy
    for i in range(trajectoryLength-slice-1, -1, -1):
        copytj[(trajectoryLength-slice)-i] = tj[i]

    # !!!
    # WARNING: TIME REVERSAL MISSING
    #can/should i reverse time?
    sim.system = copytj[-1]
    sim.system.temperature = sim.thermostat_temperature
    for j in range(slice, trajectoryLength, 1):
        sim.run()
        copytj[j] = sim.system
    return copytj

def shiftBackward(sim, tj, slice):
    """
    Perform a forward shifting move: delete a piece of trajectory
    before slice, keep the order of the trajectory and continue the
    trajectory for the missing bit on the other end.
    """
    log.debug('tps shift backward')
    trajectoryLength = len(tj)
    copytj=TrajectoryRam()
    # forward copy
    for i in range(slice, trajectoryLength, 1):
        copytj[i-slice] = tj[i]

    # continue from the end
    sim.system = copytj[-1]
    # TODO: we should reassign from the thermostat temperature like
    # sim.system.set_temperature(sim.system.thermostat.temperature)
    sim.system.temperature = sim.thermostat_temperature
    for j in range(slice):
        sim.run()
        copytj[j+(trajectoryLength-slice)] = sim.system
    return copytj


def generate_trial(sim, tj, ratio):
    """
    Generate a new trial trajectory:
    - select a random number to choose between forward and backward
      moves
    - decide then whether to shoot or to shift the trajectory
    - return the trial trajectory
    """
    r=np.random.uniform(0,1)
    if r<0.5:
    #forward
        shoot = 2.*r < ratio
        if shoot:
            trial=shootForward(sim, tj, secondHalf(tj))
        else:
            trial=shiftForward(sim, tj,firstHalf(tj))
    else:
    #backward
        shoot = 2.*(r-.5) < ratio
        if shoot:
            trial=shootBackward(sim, tj, firstHalf(tj))
        else:
            trial=shiftBackward(sim, tj, secondHalf(tj))
    return trial

def update(trajectory, attempt):
    for i in range(len(attempt)):
        trajectory[i] = attempt[i]
    # trajectory=attempt # does this work? should I also couple the simulation object?

def calculate_order_parameter(attempt):
    # extensive in time (and space)
    s = 0
    for j in range(len(attempt)):
        s += 1
    return s

def calculate_bias(tj, umbrella, k):
    """Calculate the pseudo potential, from quadratic umbrellas."""
    return 0.5*k*(calculate_order_parameter(tj)-umbrella)

def setup(bias, simulations, trajectories, umbrellas, k):
    # TODO: umbrellas should be dict of params that we unpack into calculate_bias
    for i in range(len(simulations)):
        bias[i] = calculate_bias(trajectories[i], umbrellas[i], k)

def mc_step(simulation, trajectory, biasing_field, ratio=0.25):
    """
    Perform a Monte-Carlo step in trajectory space.

    It gets in input:
    - a `simulation` object
    - a `trajectory` object related to the simulation object
    - the `umbrella` parameters
    - the `bias` array containing the value of the pseudopotential for
      the current trajectory
    - optionally, the `ratio` between shoot and shifting moves

    It returns nothing, but it should update the trajectory and the bias.
    """
    # generate trial trajectory
    q_old = calculate_order_parameter(trajectory)
    attempt = generate_trial(simulation, trajectory, ratio)
    q_new = calculate_order_parameter(attempt)
    p = np.exp(-(q_new-q_old) * biasing_field)

    if np.random.uniform() < p:
        log.debug('tps accept move bias=%s p=%s', q_new, p)
        # accept, overwrite the relevant part of trajectory
        update(trajectory,attempt)
    else:
        # reject, i.e. do nothing
        pass

    return 0.0

def mc_step_umbrella(simulation, trajectory, umbrella, k, bias, ratio=0.25):
    """
    Perform a Monte-Carlo step in trajectory space.

    It gets in input:
    - a `simulation` object
    - a `trajectory` object related to the simulation object
    - the `umbrella` parameters
    - the `bias` array containing the value of the pseudopotential for
      the current trajectory
    - optionally, the `ratio` between shoot and shifting moves

    It returns nothing, but it should update the trajectory and the bias.
    """
    # generate trial trajectory
    attempt = generate_trial(simulation, trajectory,ratio)
    attempt_bias = calculate_bias(attempt, umbrella, k)
    p = np.exp(-(attempt_bias-bias))

    if np.random.uniform()<p:
        log.debug('tps accept move bias=%s p=%s', attempt_bias, p)
        # accept, overwrite the relevant part of trajectory
        update(trajectory,attempt)
        bias = attempt_bias
    else:
        # reject, i.e. do nothing
        pass

    return bias


class TransitionPathSampling(Simulation):

    version = '%s+%s (%s)' % (__version__, __commit__, __date__)

    def __init__(self, sim, temperature, steps=0, output_path=None,
                 slices=2, k=0.01, biasing_field=0.0, restart=False):
        """
        Construct a tps instance that will run for `steps` iterations.
        
        - `sim` is either a Simulation instance or a list / tuple of
          Simulation instances.
        - `temperature` is the thermostat temperature
        - `slices` is the number of subtrajectories used to compute
        the order parameter
        - `k` is the spring constant for the umbrellas if there is
        more than one `sim` instances
        """
        Simulation.__init__(self, DryRun(), output_path=output_path,
                            steps=steps, restart=restart)
        # Non-pythonic check that sim is a list or tuple
        if not (isinstance(sim, list) or isinstance(sim, tuple)):
            self.sim = [sim]
        else:
            self.sim = sim
        # Note: the number of steps of the backend is set upon construction
        self.temperature = temperature
        self.biasing_field = [biasing_field] * len(sim)
        # Umbrellas parameters
        self.k = k  # spring constant
        self.umbrellas = range(len(self.sim))  # order parameters
        self.slices = slices

        # Trajectories objects, one per simulation instance.
        # They will have 0 frames each.
        self.trj = [TrajectoryRam() for i in range(len(self.sim))]
        # Input trajectories,
        # we only read the first frame to initialize the systems
        # TODO: it should be possible to drop this
        self.inp = [TrajectoryXYZ(self.sim[i].backend.fileinp) for i in range(len(self.sim))]
        #self.inp = [TrajectoryXYZ('data/ka_rho1.2.xyz') for i in range(len(self.sim))]

        # initial value of the bias (pseudopotential)
        # for every initial trajectory
        self.bias = np.zeros(len(self.sim))
        for sim in self.sim:
            sim.order_parameter = None

        # TODO: should we set self.backend to None??

        # Make sure base directories exist
        from atooms.core.utils import mkdir
        mkdir(self.output_path)

    def __str__(self):
        return 'Transition path sampling'

    def run_until(self, steps):
        # just shortcuts
        sim, trj = self.sim, self.trj

        # TODO: instead of run() we could use run_until() of sim[i], as we do in PT.
        # This would avoid the verbose logging... but we should use an incremental variable
        # for the running steps. Let's see how it goes...
        # Or even we could have a local logging instance as self.log which can be muted on a per simulation basis

        # TODO: FT initialise trajectories: more elegant way?
        if self.current_step == 0:
            for i in range(len(self.sim)):
                #frame 0:
                # TODO: DC this one is not necessary here, the backend has it already
                self.sim[i].system = copy(self.inp[i][0])  # copy() might not be necessary here
                # TODO: fix this hack
                self.sim[i].thermostat_temperature = self.temperature
                self.sim[i].run()
                # Trajectory frame assignement takes care of copying, so copy() is not necessary here
                self.trj[i][0] = self.sim[i].system
                # Run 0
                for j in range(1, self.slices):
                    self.sim[i].system = copy(trj[i][j-1])  # copy() might not be necessary here
                    self.sim[i].run()
                    self.trj[i][j] = self.sim[i].system
                #sim[i].order_parameter = calculate_order_parameter(trj[i])
                    
            # TPS simulation
            # calculating the value of the initial potential
            setup(self.bias, self.sim, self.trj, self.umbrellas, self.k)
            log.debug("tps bias after initialisation %s", self.bias)

        for k in range(steps - self.current_step):
            log.info('tps step %s', self.current_step + k)
            # We might have several replicas of simulations with different parameters
            for i in range(len(self.sim)): # FT: to be distributed?
                self.bias[i] = mc_step(self.sim[i], self.trj[i], self.biasing_field[i])

        # Its up to us to update our steps
        self.current_step = steps

    def _info_backend(self):
        txt = 'backend: %s' % self.sim[0]
        txt += 'output path: %s' % self.output_path
        txt += 'number of replicas: %d' % len(self.sim)
        return txt
