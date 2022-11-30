# setup.py (with automatic dependency tracking)
from setuptools import setup

setup(
    name="taffy",
    packages=["taffy", "taffy.inc", "taffy.impl", "taffy.submodules.sonLib.C.inc",
              "taffy.submodules.sonLib.C.impl", "taffy.submodules.sonLib.externalTools.cutest"],
    include_package_data=True,
    package_data={ "taffy.inc": ["*.h"],
                   "taffy.impl": ["*.c"],
                   "taffy.submodules.sonLib.C.inc": ["*.h"],
                   "taffy.submodules.sonLib.C.impl": ["*.c", "*.h"],
                   "taffy.submodules.sonLib.externalTools.cutest": ["*.h"] },
    cffi_modules=["taffy/_taffy_build.py:ffibuilder"],
)
