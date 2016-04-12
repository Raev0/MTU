EE5780 project 1-SPICE simulator in Matlab (Due on 03/17/2015)

1. You need to write your own SPICE simulator that can parse the provided netlists and solve the nonlinear DC or TR circuit analysis problems.

2. You have to use the MOSFET models described in pages 60-61, and the capacitor models described in pages 71-75 of Lecture 3.

3. The simulator output format should follow the requirements of each netlist. The following formats are used in the provided netlists:
".DC" means that you need to print out the DC solution of all varialbes.

".TR" means that you need to print out (plot) the TR solution of indicated varialbes, such as node voltage values (.PRINTNV 102),  node voltage waveforms (.PLOTNV 102), or branch currents (.PRINTBI M1). Time domain simulation specification can also be extracted from this line. For example, "TR 1.0e-11 2.0e-8" means the minimum time step size is 1.0e-11 and the end time is 2.0e-8.