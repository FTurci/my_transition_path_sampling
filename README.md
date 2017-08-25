Transition path sampling
========================

A generic, multi-core / multi-GPU implementation of the [transition path sampling](https://en.wikipedia.org/wiki/Transition_path_sampling) simulation method for atomic and molecular systems. It relies on the [atooms](https://gitlab.info-ufr.univ-montp2.fr/atooms/atooms.git) framework and on efficient simulation backends (ex. [LAMMPS](http://lammps.sandia.gov/)).

Quick start
-----------

```python
from atooms.backends.lammps import LammpsBackend
from atooms.simulation import Simulation
from atooms.transition_path_sampling import TransitionPathSampling

# Create backends and wrap them as simulation instances
file_input = 'data/lj.xyz'
cmd = """
pair_style      lj/cut 2.5
pair_coeff      1 1 1.0 1.0 2.5
neighbor        0.3 bin
neigh_modify    every 20 delay 0 check no
fix             1 all nve
"""
sim_backend = [Simulation(LammpsBackend(file_input, cmd)) for i in range(1)]
sim = [Simulation(s) for s in sim_backend]
tps = TransitionPathSampling(sim, output_path='/tmp/output_dir', steps=10)
tps.run()
```

Installation
------------
The easiest way to install sample is with pip (coming soon!)
```
pip install transition_path_sampling
```

Alternately, you can clone the code repository and install from source
```
git clone git@gitlab.info-ufr.univ-montp2.fr:daniele.coslovich/atooms/transition_path_sampling.git
cd transition_path_sampling
make install
```

Authors
-------
Daniele Coslovich: https://www.coulomb.univ-montp2.fr/perso/daniele.coslovich/

Francesco Turci: https://francescoturci.wordpress.com/
