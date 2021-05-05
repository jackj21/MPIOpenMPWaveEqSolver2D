Project 04 By: Jack Jiang, Mason Barden, Shaan Chudasama (Group 5)
==================================================================


Here is our solution to Project 4. It is divided into multiple directories
as described in the layout.


Layout
------

* The makefile and main programs are in the root directory.
* The `include` directory contains all user-written header files.
* The `src` directory contains all user-written source files.
* The `print` directory contains a Jupyter notebook and output files to read printed data.


Building
--------

To compile all programs, run

    make all

To compile all objects, run 

	make obj

To build a specific program, run (for example)

   make wave_timing


Executing
---------

Each program takes 5 command line arguments to specify the parameters of
the wave simulation. These are documented within the programs.

If you run the `wave_images` program for the required timesteps, you can run the
`plot_images.py` script to view the required plots.

To run a specific program, run (for example)

	mpirun -np <# processors> ./wave_timing N Mx My alpha nt

There are also shell script files that may be run located in the root directory.
These follow the analysis tasks accordingly by their name.

To run a shell script, run (for example)

	sbatch q1.sh
