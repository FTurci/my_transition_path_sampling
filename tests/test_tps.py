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

    def test_shift(self):
        def mobility_dist(t):
            import numpy
            from atooms.trajectory.decorators import Unfolded, filter_species
            K = []
            #for i, s in enumerate(t):
            #    print '---', i, s.particle[7].position[:], 'KKK7'
            #unfoldedtj = Unfolded(t, fixed_cm=True)
            unfoldedtj = Unfolded(t)  #, fixed_cm=True) why fixing the cm? it increases K
            # for i, s in enumerate(unfoldedtj):
            #     print '+++', i, s.particle[1].position[:]
            unfoldedtj.add_callback(filter_species, '1')
            #for i, s in enumerate(unfoldedtj):
            #     print '+++', i, s.particle[1].position[:]
            side = t[0].cell.side
            pos_0 = unfoldedtj[0].dump('pos')
            #print '+++', pos_0[7], (pos_0[7] - pos_0[7])**2, 'KKK7'
            from atooms.system.particle import cm_position
            for j in range(1, len(t)):
                pos_1 = unfoldedtj[j].dump('pos')
                #print '+++', j, pos_1[7], (pos_1[7] - pos_0[7])**2, 'KKK7'
                #if ((pos_1[k] - pos_0[k])**2>7).any():
                #    print j, k, (pos_1[k] - pos_0[k])**2
                #print 'LAMM', sorted(((pos_1 - pos_0)**2).flatten())[:3], sorted(((pos_1 - pos_0)**2).flatten())[-3:]
                K.append(numpy.sum((pos_1 - pos_0)**2))
                pos_0 = pos_1
            unfoldedtj.callbacks.pop()
            return K
        import numpy, random
        numpy.random.seed(2)
        random.seed(2)
        file_inp = 'data/ka_rho1.2.xyz'
        # Set up data structures
        cmd = """
        pair_style      lj/cut 2.5
        pair_coeff      1 1 1.0 1.0  2.5
        pair_coeff      1 2 1.5 0.8  2.0
        pair_coeff      2 2 0.5 0.88 2.2
        neighbor        0.3 bin
        neigh_modify    every 10 delay 0 check yes
        #variable t equal temp
        #variable v1 atom vx
        #fix extra all print 10 "$t" file /tmp/11
        #velocity        all create 1.0 12345
        timestep        0.004
        """
        import atooms.system
        from atooms.trajectory import TrajectoryRam
        from atooms.simulation import write_config
        trj = TrajectoryRam()
        dyn = LAMMPS(file_inp, cmd)
        dyn.system.temperature = 1.0
        dyn.system.thermostat = atooms.system.Thermostat(1.0, relaxation_time=1.0)
        block = 5000
        sim = Simulation(dyn, steps=block * 10)
        def store(sim, trj):
            trj.write(sim.system, sim.current_step)
        def check(sim):
            with open('/tmp/2.xyz', 'a') as fh:
                fh.write('{} {} {}\n'.format(sim.current_step, sim.system.temperature, sim.system.potential_energy(cache=False), [sim.system.particle[i].position[0] for i in range(2)]))
        #sim.add(check, 50)
        sim.add(store, block, trj)
        sim.temperature = 1.0  # This is needed by shift()
        sim.run()
        sim.remove(store)
        sim.remove(check)
        
        dyn = LAMMPS(file_inp, cmd)
        dyn.system.temperature = 1.0
        dyn.system.thermostat = atooms.system.Thermostat(1.0)
        sim = Simulation(dyn, steps=block)
        sim.temperature = 1.0  # This is needed by shift()
        from atooms.transition_path_sampling.core import shift_forward, shift_backward
        print numpy.mean(mobility_dist(trj))
        # new = shift_forward(sim, trj, len(trj)/2)        
        # print numpy.mean(mobility_dist(new))
        
if __name__ == '__main__':
    unittest.main(verbosity=0)

