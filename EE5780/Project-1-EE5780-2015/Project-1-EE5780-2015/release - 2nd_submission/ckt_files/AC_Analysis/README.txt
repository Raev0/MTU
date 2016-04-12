EE5780 project 1-Interconnect AC analysis (Due on 03/17/2015)


Command Line:
.AC n_point fstart fstop
        

You need to have at least one AC independent source (Vdriver1 256 0 1 in this example). The starting frequency of any AC source in the circuit is set to be fstart. The frequency increment value is defined by the increment option. Frequency is sweeped from fstart to stop. The total number of frequency points used in this simation is n_point. 


Sparse Matrices in Matlab:
For this large circuit example, you need to use sparse matrix structures to store G and C matrices. You can use matlab function spalloc(N, N, 6*N) to create an N-by-N sparse matrix with 6N initial nonzero entries.


Output:
The node voltages to be computed are shown with a ".PRINT XXX" statement. Considering ".PRINT 769", if node 256 was the input node for the circuit, then the printed output for node 769 will be the gain of your circuit from node 256 to 769. Report your amplitude and phase plots in logarithmic scale.

Model order reduction methods (AWE) can be applied for this example to obtain  extra credits.

Grading Policy:
Extra bonus credits can be obtained by completing this AC simulator.