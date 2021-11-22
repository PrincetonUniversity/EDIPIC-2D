This is a dc system periodic along the x direction, Helium plasma. 
The dc voltage (150 Volts) is applied to the top metal wall (#2), the bottom metal wall (#4) is grounded.
Both metal walls have ion induced electron emission enabled, see files init_bo_02.dat, init_bo_04.dat, and init_material_METAL2.dat .
Presently, the emission coefficient is a constant within a given energy range.
Neutral gas Helium is present, with elastic, inelastic, and ionization collisions enabled. 
Resonant charge-exchange collisions between He+ and He are included, see file init_neutral_Helium_rcx_param.dat .

The setup is configured to run on 32 MPI processes.
