#!/usr/bin/env python

"""Transition path sampling simulation driver."""

import sys
import logging

from atooms.core.utils import setup_logging
from atooms.core.utils import size, rank, barrier
from atooms.simulation.backends import LammpsBackend
from atooms.simulation import Simulation, Scheduler
from transition_path_sampling import TransitionPathSampling
from transition_path_sampling.observers import write_thermo

def main(params):
    setup_logging('atooms', level=40)
    setup_logging('transition_path_sampling', level=10)

    # TODO: in parallel, we must read input file one process at a time.
    cmd = open(params.script_file).read()
    sim = [Simulation(LammpsBackend(params.input_file, cmd),
                      steps=params.steps_backend) for i in range(params.replicas)]
    tps = TransitionPathSampling(sim, temperature=params.temperature,
                                 output_path=params.output_dir,
                                 steps=params.steps)
    tps.add(write_thermo, Scheduler(1))
    tps.restart=params.restart
    tps.run()


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', dest='verbose', action='store_true', help='verbose output')
    parser.add_argument('-T', '--temperature', dest='temperature', type=float, help='temperature')
    parser.add_argument('-n','--steps', dest='steps', type=int, default=None, help='number of TPS steps')
    parser.add_argument('-N','--steps-backend', dest='steps_backend', type=int, default=None, help='number of backend steps per TPS step')
    parser.add_argument('--replicas', dest='replicas', type=int, default=1, help='number of replicas')
    parser.add_argument('--dt', dest='dt', type=float, default=0.002, help='time step')
    parser.add_argument('--script', dest='script_file', help='script field')
    parser.add_argument('--fixcm-interval', dest='fixcm_interval', type=int, default=100, help='steps interval after which CM is fixed)')
    parser.add_argument('-i', dest='input_file', help='input_file')
    parser.add_argument('-r', dest='restart', action='store_true', help='restart')
    parser.add_argument(dest='output_dir',type=str, help='output directory')
    params = parser.parse_args()
    main(params)
