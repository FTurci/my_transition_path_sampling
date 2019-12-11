#!/usr/bin/env python

import unittest

import atooms.core
from atooms.core.utils import setup_logging
from atooms.simulation import Simulation
from atooms.backends.dryrun import DryRun
from atooms.transition_path_sampling import TransitionPathSampling

try:
    _LAMMPS = True
    from atooms.backends.lammps import LAMMPS
except ImportError:
    _LAMMPS = False

setup_logging(name='atooms.simulation', level=40)  # 20 is verbose, 40 just warnings
setup_logging(name='transition_path_sampling', level=20)
atooms.core.progress.active = False

class Test(unittest.TestCase):

    def setUp(self):
        pass

    def test_dryrun(self):
        bck = DryRun()
        sim = Simulation(bck, steps=3)
        tps = TransitionPathSampling(sim, temperature=0.8, steps=1)
        tps.run()
        self.assertEqual(tps.sim.steps, 3)
        self.assertEqual(tps.steps, 1)

    @unittest.skipIf(not _LAMMPS, 'lammps not installed')
    def test_lammps(self):
        file_inp = 'data/ka_rho1.2.xyz'
        # Set up data structures
        cmd = """
        pair_style      lj/cut 2.5
        pair_coeff      1 1 1.0 1.0  2.5
        pair_coeff      1 2 1.5 0.8  2.0
        pair_coeff      2 2 0.5 0.88 2.2
        neighbor        0.3 bin
        neigh_modify    every 20 delay 0 check no
        fix             1 all nve
        """
        sim = Simulation(LAMMPS(file_inp, cmd), steps=10)
        tps = TransitionPathSampling(sim, temperature=0.8, steps=5)
        tps.run()
        self.assertEqual(tps.sim.steps, 10)
        self.assertEqual(tps.steps, 5)

    def test_shoot_shift(self):
        frames = 10
        fileinp = 'data/ka_rho1.2.xyz'
        bck = DryRun()
        sim = Simulation(bck, steps=10)
        from atooms.transition_path_sampling import core
        from atooms.transition_path_sampling.core import shift_backward, shift_forward, \
            shoot_backward, shoot_forward, first_half, second_half
        def test_order(t):
            self.assertEqual(len(t), frames)
            return len(t)
        core.calculate_order_parameter = test_order
        tps = TransitionPathSampling(sim, temperature=0.8, steps=1, frames=frames)
        tps.run()
        nt = len(tps.trj)
        self.assertEqual(nt, len(shift_backward(tps.sim, tps.trj, first_half(tps.trj))))
        self.assertEqual(nt, len(shoot_backward(tps.sim, tps.trj, second_half(tps.trj))))
        self.assertEqual(nt, len(shift_forward(tps.sim, tps.trj, first_half(tps.trj))))
        self.assertEqual(nt, len(shoot_forward(tps.sim, tps.trj, second_half(tps.trj))))

if __name__ == '__main__':
    unittest.main(verbosity=0)

