#!/usr/bin/env python3

from setuptools import setup, find_packages
import subprocess
import distutils.command.build_py

class BuildWithMake(distutils.command.build_py.build_py):
    """
    Build using make.
    Then do the default build logic.

    """
    def run(self):
        # Call make.
        subprocess.check_call(["make"])

        # Keep building the Python stuff
        distutils.command.build_py.build_py.run(self)


setup(name="sonLib",
      version="2.0",
      description="Small general purpose library for C and Python with focus on "
      "bioinformatics.",
      author="Benedict Paten",
      author_email="benedict@soe.ucsc.edu",
      package_dir = {'': 'src'},
      packages=find_packages(where='src'))
