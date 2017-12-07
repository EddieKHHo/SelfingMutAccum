### README for simulating selfing populations experiencing fluctuating selection.

This program is written in C++. It simulates a selfing population whose genomes contain a mixture of loci experiencing time constant selection and loci experiencing selection that fluctuates over time. This program was used in the manuscript "Mutation accumulation in selfing populations under fluctuating selection" by Eddie Ho and Aneil Agrawal 2017 (unpublished).

#### BASIC OUTLINE OF SIMULATION
The program simulations diploid populations of size N assuming a life cycle of selection, reproduction and mutation. Individuals have genomes containing L loci which carry two types of loci. A proportion prc of loci experience constant selection on mutations (C-loci) and aproportion (1-pfc) of loci that experience selection that fluctuates in time (F-loci); loci are randomly placed on the genome. F-loci are further divdied into sets F1 and F2 with proportions prf0 and (1-prf0), respectively; prc + (1-prc)\*(prf0)+(1-prc)\*(1-prf0) = 1. F1 and F2 can experience different forms of fluctuating selection depending on the parameters chosen. 

Each generation, consecutive environments remained the same with probability f. With probability, 1-f, the environment was randomly selected among four possible states: {sf1(t) = s\*d0a, sf2(t) = s\*d1a }, {sf1(t) = s\*d0a, sf2(t) = s\*d1b}, {sf1(t) = s\*d0b, sf2(t) = s\*d1a}, {sf1(t) = s\*d0b, sf2(t) = s\*d1b} with probabilities (phi)\*(psi) + d, (phi)\*(1-psi) - d, (1-phi)\*(psi) - d, (1-phi)\*(1-psi) + d, respectively. Loci within F1 and F2 were expected to be under selection with strength s\*d0a and s\*d1a a proportion phi and psi of the time, respectively. Loci within F1 and F2 were expected to be under selection with strength s\*d0b and s\*d1b a proportion 1-phi and 1-psi of the time, respectively. The value of d controlls the correlation in selective states between F1 and F2 and must satisfy dmin = max {-ϕ2, -(1 - ϕ)2} when d < 0 and dmax = (1 - ϕ) ϕ when d > 0; When d = 0, selection at F1 and F2 fluctuate independently.




The correlation in selective states between set F1 and F2, ρ, was determined by the choice of d. We explore five different correlations in selective states between loci within set F1 and F2, ρ = {-1, -0.5, 0, 0.5, 1}.  When ρ = 1, this model is the same as the base model.
For simulations where selection fluctuates over time following an exponential distribution, we modeled three cases that differ in the correlation of selective states between loci within set F1 and F2. First, loci in set F1 and F2 may be completely correlated in selection (ρ = 1), therefore sf1(t) equals sf2(t) and was drawn from the exponential distribution. Second, loci in the two sets may be uncorrelated in selection (ρ = 0), therefore sf1(t) and sf2(t) were independently drawn from an exponential distribution. Lastly, loci in the two sets may be negatively correlated in selection. sf1(t) is drawn from an exponential distribution, with a value falling on the kth quantile, and sf2(t) takes on the value that is the (1-k)th quantile of the exponential distribution; this represents the minimum correlation possible between two exponentially distributed random variables (ρ = 1-π2/6 ≅ -0.645; Bladt and Neilsen 2010). 



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
* u is the diploid genome wide rate of deleterious mutation (Integer).
* h is the dominance coefficient for deleterious mutations (Integer).
* s is the unmodulated selection coefficient for deleterious mutations. The fitness effects of deleterious mutations occuring at loci under constant selection and loci under fluctuating selection are modulated by the parameters d0a, d0b, d1a, d1b and d2c. Specifically, at a locus experiencing constant selection a deleterious mutation in the heterozygous state reduces fitness by a factor of (1-h\*s\*d2c). At a locus experiencing fluctuating selection a deleterious mutation in the heterozygous state reduces fitness by a factor (1-h\*s\*d0a) and (1-h\*s\*d0a) in environment A and B, respectively. There are two sets of loci under fluctuating selection in the genome {d0a, d0b} modulates s for set F1 and {d1a, d1b} modulates s for set F2. For the simulations in the manuscript I kept d0a=d1a and d0b=d1b to make fluctuating selection equivalent between the two sets. (Float)




