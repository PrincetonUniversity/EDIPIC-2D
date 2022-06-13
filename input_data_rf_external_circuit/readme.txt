This is an example of a RF discharge with an external circuit. The simulation domain is a square, 241x241 cells. 
At the external boundary of the domain, there is a plane electrode (boundary object #1) connected in a series to the RF voltage source,
a capacitor, and the ground. There is also a C-shaped electrode (boundary obecjt #3) connected to the ground. 
Boundary segments between the two electrodes are gaps absorbing particles, 
the electrostatic potential across these gaps has a linear profile.
Parameters of the circuit, such as the rf voltage amplitude, frequency, and the capacitance of the capacitor, 
as well as the number of the electrode connected to the circuit, are defined in file init_ext_circuit.dat. 

This version of the code (13-June-2022) can simulate the circuit described above (1), as well as two other simple circuits, namely: 
(2) one metal electrode with floating potential, and (3) two metal electrodes with floating potential 
(in cases 2,3 the potentials of other electrodes must be given). To enable options 2 or 3 one must edit file pic2d_ExternalCircuits.f90 
(change 1 for 2 or 3 in line 1768). 

The methodology for the external circuit is the same as described in [Vahedi and DiPeso, JCP, 1997].  
The electric field is a superposition of electric fields created by plasma charges, electrodes with given potential, and electrodes 
connected to the circuit whose potentials are found as a solution of a system of linear equations. This system has as many equations 
as the number of electrodes connected to the circuit. The equations are defined by the circuit configuration. 
The RF discharge in this example requires only one equation for the potential of the plane RF electrode.

The circuit equation system requires the charge of each electrode connected to the circuit to be represented as a sum of terms 
proportional to the potentials of these electrodes. The code has all the necessary procedures and data structures in place, 
see file pic2d_ExternalCircuits.f90. Producing the system of equations for a specific circuit is the responsibility of the user,
see subroutine SOLVE_EXTERNAL_CONTOUR in pic2d_ExternalCircuits.f90.

Time dependences of circuit parameters, such as the electrode potential, charge, the circuit current, and some others are saved in 
history_ext_circuit.dat, see lines 1863-1875 in file pic2d_ExternalCircuits.f90. If one uses a more complex circuit, this 
diagnostics may have to be updated as well.

If the code does not find file init_ext_circuit.dat, simulation continues without the external circuit. 
If the simulation has to be carried out without thet circuit, make sure that file init_ext_circuit.dat 
is not in the simulation directory.

Note that format of input data file init_ext_circuit.dat has changed compared to the older (March 2022) version.  
Parameters of the external circuit power supply must be specified in a single line. It is possible to introduce a power supply
with a constant voltage, a harmonically oscillating voltage, a voltage with the time profile defined by a waveform data file, 
or a combination of the above. The amplitude of the oscillating part may be a separate function of time defined via 
an amplitude profile data file. If the code does not find the waveform and/or the amplitude profile data file for the 
corresponding external circuit power supply, the voltage of the power supply will be only as specified in the 
init_ext_circuit.dat (constant + sinusoidal oscillations).

Note that while the new input data file init_ext_circuit.dat allows to introduce numerous power supplies,
resistors, capacitors, and inductors in the external circuit, the code still can solve only the circuit with one HF power source 
and one capacitor as described above. For a different circuit, the system of equations solved in subroutine 
SOLVE_EXTERNAL_CONTOUR in pic2d_ExternalCircuits.f90 must be modified.

An example of the external circuit with the waveforms and the amplitude profiling is in
https://github.com/PrincetonUniversity/EDIPIC-2D/tree/main/input_data_rf_external_circuit_waveform_pulsed_amplitude .

Format of the waveform and amplitude profile data files is similar to that of the corresponding files 
for the potential of inner/boundary objects.

The example is configured to run with 16 MPI processes.

