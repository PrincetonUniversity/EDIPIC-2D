This is a plasma periodic along the x-direction bounded between two metal electrodes in the y-direction. 
One electrode (boundary object #3) is grounded.
The electrostatic potential of the other electrode (#1) is a prescribed function of time. 
Actually, it is a product of two functions.
One function gives high frequency oscillations with the shape defined in file init_bo_01_waveform.dat.
The other function controls variation of the amplitude of the oscillating potential with time. 
The shape of this function is defined in file init_bo_01_amplitude_profile.dat. The amplitude profile function is also periodic.
Its period is the time of the last point of the shape function defined in  init_bo_01_amplitude_profile.dat.

If the code does not find file init_bo_NN_amplitude_profile.dat, 
then the amplitude of the oscillating potential of boundary/inner object #NN (if requested) is constant in time.

The simulation runs on 16 MPI processes.
