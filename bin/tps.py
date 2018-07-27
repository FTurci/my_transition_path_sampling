#!/usr/bin/env python

"""Transition path sampling simulation driver."""

if __name__ == '__main__':
    from atooms.transition_path_sampling.api import main
    from argh import dispatch_command
    dispatch_command(main)

