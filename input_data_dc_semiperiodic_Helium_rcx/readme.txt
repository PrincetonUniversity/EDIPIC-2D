This is a dc system periodic along the x direction, Helium plasma. 
The dc voltage (150 Volts) is applied to the top metal wall (#2), the bottom metal wall (#4) is grounded.
Both metal walls have ion induced electron emission enabled, see files init_bo_02.dat, init_bo_04.dat, and init_material_METAL2.dat .
Presently, the emission coefficient is a constant within a given energy range.
Neutral gas Helium is present, with elastic, inelastic, and ionization collisions enabled. 
Resonant charge-exchange collisions between He+ and He are included, see file init_neutral_Helium_rcx_param.dat .

The vacuum permittivity eps_0 is acquired from file init_physconstants.dat, in this test case it is equal to the true value (8.85e-12 F/m).
Note, if file init_physconstants.dat is not found, by default the code sets the value of eps_0_Fm equal to the true value.

Also, the simulation will create checkpoints using ordinary Fortran I/O commands (below referred to as POSIX) instead of MPI I/O, with multiple binary files (one per MPI process) for each checkpoint.
The checkpoint files are placed into directories with names checkdir_TTTTTTTT where TTTTTTTT is the eight digit time step of the checkpoint.
To switch between creating MPIIO checkpoints (each checkpoint is one large binary file) and POSIX checkpoints, use positive/negative values of  
the number of ion cycles between checkpoints in file init_simcontrol.dat, respectively (presently it is -100).

The process of restarting does not change - just specify the timestep of the checkpoint to be used and make sure that either MPIIO or POSIX checkpoint with this number exists.

The setup is configured to run on 32 MPI processes.
