This is a dc plasma discharge. The simulation domain is a rectangle. The external boundary is grounded metal. 
In the center of the domain there is a small square metal electrode with the potential of -200 Volts. 
The surface of the inner electrode emits secondary electrons if it is bombarded by energetic ions.

The diagnostics output includes snapshots averaged over multiple time steps 
(averaging time interval is about 10 ns in this case). New parameters saved in snapshots are 
electron and ion heat flow vector components
defined as q_s^vector = 1/2 n_s m_s <c_s^2 c_s^vector>, where c_s^vector = v_s^vector - u_s^vector 
is the random velocity, u_s^vector = <v_s^vector> is the drift velocity, <> is averaging over particles, 
subscript s denotes species, see [Schunk, Reviews of Geophysics and Sppace Physics, vol.15, p.429, 1977].

Control of creation of time-averaged snapshots is via input file init_avgsnapshots.dat.

If file init_avgsnapshots.dat is not found, the simulation continues without creating time-averaged snapshots.

Note that one can request instant snapshots of the heat flow vector components 
(in file init_snapshots.dat). Most likely, however, such a snapshots will be very noisy, which is why
it is recommended to use the time averaged snapshots.

The simulation is configured for 16 MPI processes.
