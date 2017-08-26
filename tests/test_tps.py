#!/usr/bin/env python

import unittest

from atooms.utils import setup_logging
from atooms.simulation import Simulation
from atooms.simulation.backend_dryrun import DryRunBackend
from transition_path_sampling import TransitionPathSampling

try:
    _LAMMPS = True
    from atooms.simulation.backend_lammps import LammpsBackend
except ImportError:
    _LAMMPS = False

setup_logging(name='atooms.simulation', level=40)  # 20 is verbose, 40 just warnings
setup_logging(name='transition_path_sampling', level=20)

class Test(unittest.TestCase):

    def setUp(self):
        pass

    def test_dryrun(self):
        n = 1
        fileinp = 'data/ka_rho1.2.xyz'
        bck = [DryRunBackend() for i in range(n)]
        # TODO: temporarily needed by tps, should be removed
        for b in bck:
            b.fileinp = fileinp
        sim = [Simulation(bck[i], steps=3) for i in range(n)]
        tps = TransitionPathSampling(sim, temperature=0.8, steps=1)
        tps.run()
        self.assertEqual(tps.sim[0].steps, 3)
        self.assertEqual(tps.steps, 1)

    @unittest.skipIf(not _LAMMPS, 'lammps not installed')
    def test_lammps(self):
        n = 1
        file_inp = 'data/ka_rho1.2.xyz'
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
        sim = [Simulation(LammpsBackend(file_inp, cmd), steps=10) for i in range(n)]
        tps = TransitionPathSampling(sim, temperature=0.8, steps=5)
        tps.run()
        self.assertEqual(tps.sim[0].steps, 10)
        self.assertEqual(tps.steps, 5)

if __name__ == '__main__':
    unittest.main(verbosity=0)

