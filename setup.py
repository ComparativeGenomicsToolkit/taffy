# setup.py (with automatic dependency tracking)
from setuptools import setup

setup(
    name="taffy",
    packages=["taffy", "taffy.inc", "taffy.impl", "taffy.subModules.sonLib.C.inc",
              "taffy.subModules.sonLib.C.impl", "taffy.subModules.sonLib.externalTools.cutest"],
    include_package_data=True,
    cffi_modules=["taffy/_taffy_build.py:ffibuilder"],
)
