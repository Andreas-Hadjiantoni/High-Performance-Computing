# High Performance Computing

The task of this coursework was to optimize some stencil operations on an image.
The optimizations were first done on a single core. Then, multiple cores were incorporated.

## Single core optimizations:

- Vectorisation of operations.
- Optimised the ordering of operations for better cache exploitation.
- Optimised code so that the number of "Branch" instructions in the main loop of the executable is minimized. 

## Multiple Cores:

- Message Passing Interface(MPI) was used to employ multiple worker cores. Each core was given a region in the image, and the necessary info was communicated between cores.
- Three main modes/architectures were implemented:
    1) Synchronous data exchnge.
    2) Asynchronous data exchnge.
    3) Synchronous data Exchange with cartesian Topology.

Intel Advisor was used throughout the development for code profiling. The code was run on BlueCrystal, the University Of Bristol's supercomputer:

https://www.acrc.bris.ac.uk/acrc/phase3.htm

This image is an example input for the program:

![](stencilInput.pgm)

And this is the output:

![](stencilOutput.pgm)
