This is an electron beam-plasma system. The system is periodic along the x direction (the periodic boundaries are boundary objects 1 and 3). 
The top boundary is a metal wall (boundary object 2) which does not emit electrons.
The bottom boundary is split into three parts. Two parts (boundary objects 4 and 6) are metal walls without emission.
One part (boundary object 5) is a metal wall, also, an electron beam is injected through this boundary.
Boundary object 5 is between boundary objects 4 and 6. All metal walls have zero potential.
Since this is a system periodic in one direction and no inner objects are included, the FFT-based field solver is used.
Neutral gas Argon is present, with elastic, excitation, and ionization collisions turned on.

The setup is configured to run with 32 MPI processes.

