#!/usr/bin/env python

import os

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

with open('README.md') as f:
    readme = f.read()

with open('atooms/transition_path_sampling/_version.py') as f:
    exec(f.read())

setup(name='transition_path_sampling',
      version=__version__,
      description='An atooms package for transition path sampling simulations',
      long_description=readme,
      author='Daniele Coslovich',
      author_email='daniele.coslovich@umontpellier.fr',
      url='https://gitlab.info-ufr.univ-montp2.fr/daniele.coslovich/atooms/transition_path_sampling',
      packages=['transition_path_sampling'],
      license='GPLv3',
      install_requires=['atooms'],
      scripts=['bin/tps.py'],
      classifiers=[
          'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
          'Programming Language :: Python',
      ]
     )
