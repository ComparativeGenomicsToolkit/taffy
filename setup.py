# setup.py (with automatic dependency tracking)
from setuptools import setup

setup(
    name="taffy",
    packages=["taffy", "taffy.inc", "taffy.impl", "taffy.subModules.sonLib.C.inc",
              "taffy.subModules.sonLib.C.impl", "taffy.subModules.sonLib.externalTools.cutest"],
    include_package_data=True,
    setup_requires=["cffi>=1.0.0"],
    cffi_modules=["taffy/_taffy_build.py:ffibuilder"],
    install_requires=["cffi>=1.0.0"],
)
