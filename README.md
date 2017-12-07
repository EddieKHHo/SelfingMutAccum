### README for simulating selfing populations experiencing fluctuating selection.

This program is written in C++. It simulates a selfing population whose genomes contain a mixture of loci experiencing time constant selection and loci experiencing selection that fluctuates over time. This program was used in the manuscript "Mutation accumulation in selfing populations under fluctuating selection" by Eddie Ho and Aneil Agrawal 2017 (unpublished).

#### RUNNING THE PROGRAM  
This program is meant to be run on a command line. 
It must first be compiled with these files all within the same folder and typing './make':
* GlobalConstants.h
* IncludeFiles.h
* SharedFunctions.h
* SharedFunction.cpp
* InputPar.h
* InputPar.cpp
* Indiv.h
* Indiv.cpp
* Popul.h
* Popul.cpp
* World.h
* World.cpp
* main.cpp
* MakeFile

After compiling the program is run by reading a text file containing all the necessary parameters following this syntax:

Below is an example of a parameter text file (this is also provided as SampleParameters.txt):  
-dir /User/OutputResult/  
-mfn Simulation1  
-nr 10  
-mng 450000  
-xmin 10000  
-xdat 100 
-qneu 1  
-l 5000  
-n 10000  
-u 0.5  
-h 0.25  
-s 0.05  
-ub 0.001  
-hb 0.5  
-sb 0.1  
-is 0.98  
-ir 1  
-sto 1  
-f 0.95  
-phi 0.5  
-psi 0.5  
-eps 0  
-eta 0  
-d0a 1  
-d0b 0  
-d1a 1  
-d1b 0  
-d2c 0.5  
-prc 1  
-prf0 0.5  
#No blank lines for last line  

Definition of each parameter:
* dir is the directory for outputing the result files for the simulation.
* mfn is the file name for the simulation being run (Integer).
* nr is the number of replicates to run the simulation (Integer).
* mng is the number of generation to run the simulation before it ends (Integer).
* xmin is the generation from which data collection will start (Integer).
* xdat is the interval, in generations, between collection data (Integer).
* qneu indicates whether the simulation should track neutral loci for measuring effective populations size (0 = No, 1 = Yes).
* l is the number of loci in the genome (Integer).
* n is the population size of the selfing population (Integer).
* u is the diploid genome wide rate of mutation (Integer).
* h is the dominance coefficient for deleterious mutations (Integer).
* s is the unmodulated selection coefficient for deleterious mutations. The effects of these deleterious mutations occuring at loci under constant selection and loci under fluctuating selection is modulated by the parameters d0a, d0b, d1a, d1b and d2c. Specifically at a locus experiencing constant selection a deleterious mutation in the heterozygous state reduces fitness by a factor of (1-h*s*d2c). At a locus experiencing fluctuating selection a deleterious mutation in the heterozygous state reduces fitness by a factor (1-h*s*d0a) and (1-h*s*d0a) in environment A and B, respectively. There are two sets of loci under fluctuating selection in the genome {d0a, d0b} modulates s for set F1 and {d1a, d1b} modulates s for set F2. For the simulations in the manuscript I kept d0a=d1a and d0b=d1b to make fluctuating selection at the two sets equivalent.




