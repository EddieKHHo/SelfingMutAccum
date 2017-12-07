README for simulating selfing populations experiencing fluctuaing selection.

This program is written in C++. It simulates a selfing population whose genomes contain a mixture of loci experiencing time constant selection and loci experiencing selection that fluctuates over time. This program was used in the manuscript "Mutation accumulation in selfing populations under fluctuating selection" by Eddie Ho and Aneil Agrawal 2017 (unpublished).

RUNNING THE PROGRAM
This program is meant to be run on a command line. It must first be compiled with these files all within the same folder and typing './make':
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
-xdat 100  
-xmin 10000  
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

