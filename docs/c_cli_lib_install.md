# Installing Taffy CLI/C Library

Do build this repo clone the repo as follows and then make:

    git clone https://github.com/ComparativeGenomicsToolkit/taffy.git --recursive
    cd taffy && make

To test the installation do:

    make test

This will run the unit tests. You should see that all tests pass okay. You will
then want to add the taf/bin directory to your path. 