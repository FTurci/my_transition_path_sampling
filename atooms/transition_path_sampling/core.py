"""Sketch of TPS simulation."""

import logging
import numpy as np
from copy import copy
from atooms.trajectory import TrajectoryRam, TrajectoryXYZ
from atooms.simulation import Simulation
from atooms.backends.dryrun import DryRun
from atooms.transition_path_sampling import __version__, __date__, __commit__

log = logging.getLogger(__name__)


def first_half(tj, margin=0):
    """Return an integer in [margin,last/2)"""
    last = len(tj)
    return margin+np.random.randint(last/2-margin)

def second_half(tj, margin=0):
    """Return an integer in [last/2,last-margin)"""
    last = len(tj)
    return last-margin-1-np.random.randint(last/2-margin)

def shoot_forward(sim, tj, frame):
    """
    Perform a _forward shooting move: from a selected frame, randomise
    the velocities and integrate a new piece of trajectory _forward in
    time.
    """
    log.debug('tps shoot forward')
    sim.system = tj[frame]
    # We completely reshuffle the velocities, which is a big perturbation...
    sim.system.set_temperature(sim.temperature)

    # overwrtiting from frame till the end
    for j in range(frame, len(tj)):
        sim.run()
        tj[j] = sim.system
    return tj

def shoot_backward(sim, tj, frame):
    """
    Perform a backward shooting move: from a selected frame, randomise
    the velocities and integrate a new piece of trajectory backward in
    time.
    """
    log.debug('tps shoot backward')
    sim.system = tj[frame]
    sim.system.set_temperature(sim.temperature)
    # !!!
    # Possible issue with time reversal
    
    for j in range(frame-1, -1, -1):
        sim.run()
        tj[j] = sim.system
    return tj

def shift_forward(sim, tj, frame):
    """
    Perform a _forward shifting move: delete a piece of trajectory
    starting from frame reverse the order of the trajectory and
    continue the trajectory for the missing bit on the other end.
    """
    log.debug('tps shift _forward')
    last = len(tj)
    copytj = TrajectoryRam()
    # reverse copy
    for i in range(last-frame-1, -1, -1):
        # copytj[(last-frame)-i] = tj[i]
        system = tj[i]
        for p in system.particle:
            p.velocity *= -1
        copytj[(last-frame)-i] = system

    # !!!
    # Possible issue with time reversal
    #sim.system.set_temperature(sim.temperature)
    sim.system = copytj[-1]
    sim.system.set_temperature(sim.temperature)

    for j in range(frame, last, 1):
        sim.run()
        copytj[j] = sim.system
    return copytj

def shift_backward(sim, tj, frame):
    """
    Perform a _forward shifting move: delete a piece of trajectory
    before frame, keep the order of the trajectory and continue the
    trajectory for the missing bit on the other end.
    """
    log.debug('tps shift backward')
    last = len(tj)
    copytj=TrajectoryRam()
    # _forward copy
    for i in range(frame, last, 1):
        copytj[i-frame] = tj[i]

    # continue from the end
    sim.system = copytj[-1]
    #
    sim.system.set_temperature(sim.temperature)
    for j in range(frame):
        sim.run()
        copytj[j+(last-frame)] = sim.system
    return copytj


def generate_trial(sim, tj, ratio):
    """
    Generate a new trial trajectory:
    - select a random number to choose between _forward and backward
      moves
    - decide then whether to shoot or to shift the trajectory
    - return the trial trajectory
    """
    r=np.random.uniform(0,1)
    if r<0.5:
    #_forward
        shoot = 2.*r < ratio
        if shoot:
            trial=shoot_forward(sim, tj, second_half(tj))
        else:
            trial=shift_forward(sim, tj,first_half(tj))
    else:
    #backward
        shoot = 2.*(r-.5) < ratio
        if shoot:
            trial=shoot_backward(sim, tj, first_half(tj))
        else:
            trial=shift_backward(sim, tj, second_half(tj))
    return trial

def update(trajectory, attempt):
    for i in range(len(attempt)):
        trajectory[i] = attempt[i]

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

def mc_step(simulation, trajectory, biasing_field,ratio=0.25):
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
    q_old = simulation.order_parameter 

    attempt = generate_trial(simulation, trajectory, ratio)
    q_new = calculate_order_parameter(attempt)
    p = np.exp(-(q_new-q_old) * biasing_field)

    if np.random.uniform() < p:
        log.debug('tps accept move bias=%s p=%s', q_new, p)
        # accept, overwrite the relevant part of trajectory
        update(trajectory,attempt)
        simulation.order_parameter = q_new
        return 1.0

    else:
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
        
        return 1.0
    else:
        # reject, i.e. do nothing
        return 0.0

class TransitionPathSampling(Simulation):

    version = '%s+%s (%s)' % (__version__, __commit__, __date__)

    def __init__(self, sim, temperature, steps=0, output_path=None,
                 frames=2, k=0.01, biasing_field=0.0, restart=False,
                 shift_weight=1, shoot_weight=1):
        """
        Construct a tps instance that will run for `steps` iterations.
        
        - `sim` is either a Simulation instance or a list / tuple of
          Simulation instances.
        - `temperature` is the thermostat temperature
        - `frames` is the number of subtrajectories used to compute
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
        self.biasing_field = [biasing_field] * len(self.sim)
        # Umbrellas parameters
        self.k = k  # spring constant
        self.umbrellas = range(len(self.sim))  # order parameters
        self.frames = frames
        self.acceptance = 0.0
        # Ratio of moves
        self._ratio = [shoot_weight / float(shoot_weight + shift_weight),
                       shift_weight / float(shoot_weight + shift_weight)]

        # Internal trajectory objects, one per simulation instance
        self.trj = [TrajectoryRam() for i in range(len(self.sim))]

        # Initial value of the bias (pseudopotential)
        # for every initial trajectory
        self.bias = np.zeros(len(self.sim))
        for sim in self.sim:
            sim.order_parameter = None
            sim.temperature = self.temperature

        # TODO: might be needed for PT
        # Make sure base directories exist
        # from atooms.core.utils import mkdir
        # mkdir(self.output_path)

    def __str__(self):
        return 'Transition path sampling'

    def run_until(self, steps):
        # Just shortcuts
        sim, trj = self.sim, self.trj

        # Initialize trajectories
        if self.current_step == 0:
            log.debug("generating first trajectory")
            for i in range(len(self.sim)):
                # Trajectory frame assignement takes care of copying
                self.trj[i][0] = self.sim[i].system
                # Run 0
                for j in range(1, self.frames):
                    # self.sim[i].system = copy(trj[i][j-1])  # copy() might not be necessary here
                    self.sim[i].run()
                    self.trj[i][j] = self.sim[i].system

                q_old = calculate_order_parameter(self.trj[i])
                sim[i].order_parameter = q_old
                    
            # TPS simulation
            # Calculate the value of the initial potential
            setup(self.bias, self.sim, self.trj, self.umbrellas, self.k)
            log.debug("tps bias after initialisation %s", self.bias)

        for k in range(steps - self.current_step):
            log.debug('tps step %s', self.current_step + k)           
            # We might have several replicas of simulations with different parameters
            for i in range(len(self.sim)): # FT: to be distributed?
                self.sim[i].system.set_temperature(self.temperature)
                self.acceptance += mc_step(self.sim[i], self.trj[i],
                                           self.biasing_field[i], ratio=self._ratio[0])

        # It's up to us to update our steps
        self.current_step = steps

    def _info_backend(self):
        txt = 'backend: %s' % self.sim[0]
        txt += 'output path: %s' % self.output_path
        txt += 'number of replicas: %d' % len(self.sim)
        return txt
