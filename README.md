# High Performance Computing

The task of this coursework was to optimize some stencil operations on a 2D grid.
The optimizations were first done on a single core. Then, multiple cores were incorporated.

## Single core optimizations:

- Vectorisation of operations
- Optimised the ordering of operations for better cache exploitation

## Multiple Cores:

- Message Passing Interface(MPI) was used to employ multiple worker cores. Each core was given a region in the grid, and the necessary info was communicated between cores.
- Three main modes/architectures were implemented:
    1) Synchronous data exchnge
    2) Asynchronous data exchnge
    3) Synchronous data Exchange with cartesian Topology

The code was run on BlueCrystal, the University Of Bristol's supercomputer:

https://www.acrc.bris.ac.uk/acrc/phase3.htm
