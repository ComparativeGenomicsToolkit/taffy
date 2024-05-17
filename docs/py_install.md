# Installing Python Library

Tl;dr:

```
pip install taffy
```

Longer version: to avoid problems with conflicting versions of dependencies on your system, we strongly recommend installing
the package inside a Python 3 [virtual environment](https://virtualenv.pypa.io/en/stable/).
To install the `virtualenv` command, if you don't have it already, run:

```
python3 -m pip install virtualenv
```

To set up a virtual environment in the directory `taffy_env`, run:

```
python3 -m virtualenv -p python3.XX taffy_env
```

Where XX is the specific version of Python3 that you wish to use (you can omit the .XX if you want to use the default). Also, note that I have tested this with 3.9.

Then, to enter the virtualenv, run:

```
source taffy_env/bin/activate
```

You can always exit out of the virtualenv by running `deactivate`.
Finally, install taffy with pip:

To install these notebooks in Python, clone the repo:
```
pip install taffy
```

To check that it worked try launching a python interpretor and importing
taffy.lib.

# Installing Python Library From Source

To build the Python library you must first install [htslib](http://www.htslib.org/) for bgzip support.

If you are building the C library from source you can also build the Python library as follows:

```
git clone https://github.com/benedictpaten/taf.git --recursive
cd taf && make test
```

This will build and the test the C installation. Then create a virtualenv
if you haven't already:

```
python3 -m pip install virtualenv
python3 -m virtualenv -p python3.XX taffy_env
source taffy_env/bin/activate
```

You can then build from source and test the distribution by running:

```
python3 -m pip install build
python3 -m build
pip install .
cd tests && python3 taffyTest.py
```

All these commands should succeed. Note, because of the current package structure
trying to import the library from the root directory will fail.