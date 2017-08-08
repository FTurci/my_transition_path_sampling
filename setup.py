#!/usr/bin/env python

import os

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

with open('README.md') as f:
    readme = f.read()

setup(name='sample',
      version='0.1',
      description='',
      long_description=readme,
      author='Daniele Coslovich',
      author_email='daniele.coslovich@umontpellier.fr',
      url='https://gitlab.info-ufr.univ-montp2.fr/daniele.coslovich/sample',
      packages=[],
      #packages=find_packages(exclude=('tests', 'docs')),
      license='GPLv3',
      install_requires=[],
      scripts=[],
      classifiers=[
          'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
          'Programming Language :: Python',
      ]
     )
