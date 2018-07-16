# This file is part of atooms
# Copyright 2010-2014, Daniele Coslovich

import logging
logging.getLogger(__name__).addHandler(logging.NullHandler())

try:
    from _version import __version__
except ImportError:
    __version__ = ""

try:
    from _commit import __commit__, __date__
except ImportError:
    __commit__ = ""
    __date__ = ""

from .core import TransitionPathSampling
