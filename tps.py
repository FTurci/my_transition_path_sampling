import logging
import numpy as np
from atooms.simulation import Simulation
from atooms.trajectory import TrajectoryRam

_log = logging.getLogger(__name__)

#all this should be organised in a class
runsteps = 10
margin = 0
temperature = 1.0
#spring constant for the umbrellas
k = 0.01

def firstHalf(tj):
    """Return an integer in [margin,trajectoryLength/2)"""
    trajectoryLength = len(tj)
    return margin+np.random.randint(trajectoryLength/2-margin)

def secondHalf(tj):
    """Return an integer in [trajectoryLength/2,trajectoryLength-margin)"""
    trajectoryLength = len(tj)
    return trajectoryLength-margin-1-np.random.randint(trajectoryLength/2-margin)

def shootForward(sim, tj, slice):
    """
    Perform a forward shooting move: from a selected slice, randomise
    the velocities and integrate a new piece of trajectory forward in
    time.
    """
    _log.debug('tps shoot forward')
    sim.system = tj[slice]
    sim.system.temperature = temperature

    # overwrtiting from slice till the end
    for j in range(slice, len(tj)):
        sim.run(runsteps)
        tj[j] = sim.system
    return tj

def shootBackward(sim, tj,slice):
    """
    Perform a backward shooting move: from a selected slice, randomise
    the velocities and integrate a new piece of trajectory backward in
    time.
    """
    _log.debug('tps shoot backward')
    sim.system = tj[slice]
    sim.system.temperature = temperature
    # !!!
    # WARNING: TIME REVERSAL MISSING
    #can/should i reverse time? 
    for j in range(slice-1, -1, -1):
        sim.run(runsteps)
        tj[j] = sim.system
    return tj

def shiftForward(sim, tj,slice):
    """
    Perform a forward shifting move: delete a piece of trajectory
    starting from slice reverse the order of the trajectory and
    continue the trajectory for the missing bit on the other end.
    """
    _log.debug('tps shift forward')
    trajectoryLength = len(tj)
    copytj = TrajectoryRam()
    # reverse copy
    for i in range(trajectoryLength-slice-1, -1, -1):
        copytj[(trajectoryLength-slice)-i] = tj[i]

    # !!!
    # WARNING: TIME REVERSAL MISSING
    #can/should i reverse time? 
    sim.system = copytj[-1]
    sim.system.temperature = temperature
    for j in range(slice, trajectoryLength, 1):
        sim.run(runsteps)
        copytj[j] = sim.system
    return copytj
    
def shiftBackward(sim, tj, slice):
    """
    Perform a forward shifting move: delete a piece of trajectory
    before slice, keep the order of the trajectory and continue the
    trajectory for the missing bit on the other end.
    """
    _log.debug('tps shift backward')
    trajectoryLength = len(tj)
    copytj=TrajectoryRam()
    # forward copy
    for i in range(slice, trajectoryLength, 1):
        copytj[i-slice] = tj[i]

    # continue from the end
    sim.system = copytj[-1]
    sim.system.temperature = temperature
    for j in range(slice):
        sim.run(runsteps)
        copytj[j] = sim.system
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
    trajectory=attempt # does this work? should I also couple the simulation object?

def calculate_order_parameter(attempt):
    # extensive in time (and space)
    s = 0
    for j in range(len(attempt)):
        s += 1
    return s

def calculate_bias(tj, umbrella):
    """Calculate the pseudo potential, from quadratic umbrellas."""
    return 0.5*k*(calculate_order_parameter(tj)-umbrella)

def setup(bias,simulations,trajectories,umbrellas):
    for i in range(len(simulations)):
        bias[i] = calculate_bias(trajectories[i], umbrellas[i])

def mc_step(simulation, trajectory, umbrella, bias, ratio=0.25):
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
    attempt_bias = calculate_bias(attempt, umbrella)
    p = np.exp( -(attempt_bias-bias) )

    if np.random.uniform()<p:
        _log.debug('tps accept move bias=%s p=%s', attempt_bias, p)
        # accept, overwrite the relevant part of trajectory
        update(trajectory,attempt)
        bias = attempt_bias
    else:
        # reject, i.e. do nothing
        pass
