### README for simulating selfing populations experiencing fluctuating selection.

This program is written in C++. It simulates a selfing population whose genomes contain a mixture of loci experiencing time constant selection and loci experiencing selection that fluctuates over time. This program was used in the manuscript "Mutation accumulation in selfing populations under fluctuating selection" by Eddie Ho and Aneil Agrawal 2017 (unpublished).

#### BASIC OUTLINE OF SIMULATION
The program simulations diploid populations of size *N* assuming a life cycle of selection, reproduction and mutation. Individuals have genomes containing *L* loci which carry two types of loci. A proportion prc of loci experience constant selection on mutations (C-loci) and aproportion (1-*pfc*) of loci that experience selection that fluctuates in time (F-loci); loci are randomly placed on the genome. F-loci are further divdied into sets F1 and F2 with proportions prf0 and (1-*prf0*), respectively; prc + (1-*prc*)\*(*prf0*)+(1-*prc*)\*(1-*prf0*) = 1. F1 and F2 can experience different forms of fluctuating selection depending on the parameters chosen. The fitness of an individual was calculated assuming loci act multiplicatively as *W(t)* = (1 - *s d2c*)<sup>xc</sup> (1 - *h sd2c*)<sup>yc</sup> (1 - *sf1(t)*)<sup>xf1</sup> (1 - *h sf1(t)*)<sup>yf1</sup> (1 - *s sf2(t)*)<sup>xf2</sup> (1 - *h sf2(t)*)<sup>yf2</sup>. *h* and *s* is the dominance and the unmodulated selection cofficient of deleterious mutations, respectively. xc and yc represents the number C-loci that are homozygous and heterozygous mutants, xf1 and yf1 represents the number of F1-loci are homozygous and heterozygous mutants and xf2 and yf2 represents the number of F1-loci are homozygous and heterozygous mutants. If beneficial mutations exist at C-loci then fitness calcualted with an additional multiplicative factor of (1 - *sb*)<sup>xb</sup> (1 - *hb sb*)<sup>yb</sup>. *hb* and *sb* is the dominance and the selection cofficient of beneficial leterious mutations, respectively. xb and yb are the numbers of C-loci that are homozygous and heterozygous for the beneficial mutations, respectively. 

Each generation, consecutive environments remained the same with probability *f*. With probability, 1-*f*, the environment was randomly selected among four possible states: {*sf1(t)* = *s d0a*, *sf2(t)* = *s d1a* }, {*sf1(t)* = *s d0a*, *sf2(t)* = *s d1b*}, {*sf1(t)* = *s d0b*, *sf2(t)* = *s d1a*}, {*sf1(t)* = *s d0b*, *sf2(t)* = *s d1b*} with probabilities (*phi*)(*psi*) + *eps*, (*phi*)(1-*psi*) - *eps*, (1-*phi*)(*psi*) - *eps*, (1-*phi*)(1-*psi*) + *eps*, respectively. Loci within F1 and F2 were expected to be under selection with strength s\*d0a and s\*d1a a proportion *phi* and *psi* of the time, respectively. Loci within F1 and F2 were expected to be under selection with strength *s d0b* and *s d1b* a proportion 1-*phi* and 1-*psi* of the time, respectively. The value of *eps* controls the correlation in selective states between F1 and F2 and must satisfy eps_min = max {-(*phi*)(*psi*), -(1-*phi*)(1-*psi*)} when *eps* < 0 and eps_max = min{(*phi*)(1-*psi*), (1-*phi*)(*psi*)} when *eps* > 0; when *eps* = 0, selection at F1 and F2 fluctuate independently.

The parameter sto controls how selection fluctuates stochastically. When sto = 1, selection fluctuates as described and only takes on two values for each set of F-loci (i.e. sf1(t) = {*s d0a*, *s d0b*}, sf2(t) = {*s d1a*, *s d1b*}; I refer to this as the bimodal model. When sto = 2, the selection coefficient is sampled from an exponentail distribution with a mean value of s when fluctuations occurs; I refer to this as the exponential model. We allowed for three cases that differ in the correlation of selective states between loci within set F1 and F2 under the exponential model. First, loci in set F1 and F2 may be completely correlated in selection (*eta* = 1), therefore *sf1(t)* equals *sf2(t)* and was drawn from the exponential distribution. Second, loci in the two sets may be uncorrelated in selection (*eta* = 0), therefore *sf1(t)* and *sf2(t)* were independently drawn from an exponential distribution. Lastly, loci in the two sets may be negatively correlated in selection (*eta* = -1). *sf1(t)* is drawn from an exponential distribution, with a value falling on the kth quantile, and *sf2(t)* takes on the value that is the (1-k)th quantile of the exponential distribution; this represents the minimum correlation possible between two exponentially distributed random variables (correlation = 1-(π<sup>2</sup>)/6 ≅ -0.645; Bladt and Neilsen 2010). 

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
./main /ParameterFolder/SampleParameters.txt

#### EXAMPLE PARAMETER FILE
Below is an example of a parameter text file (this is also provided as SampleParameters.txt):  
-dir /User/OutputResult/  
-mfn MySim  
-nr 10  
-mng 60000  
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
-d2c 0.506  
-prc 1  
-prf0 0.5  

Definition of parameters:
* dir is the directory for outputing the result files for the simulation.
* mfn is the file name for the simulation being run (Integer).
* nr is the number of replicates to run the simulation (Integer).
* mng is the number of generation to run the simulation before it ends (Integer).
* xmin is the generation from which data collection will start (Integer).
* xdat is the interval, in generations, between collection data (Integer).
* qneu indicates whether the simulation should track neutral loci for measuring effective populations size (0 = No, 1 = Yes).
* l is the number of loci in the genome (Integer).
* n is the population size of the selfing population (Integer).
* u is the diploid genome wide rate of deleterious mutation (Float).
* h is the dominance coefficient for deleterious mutations (0 to 1).
* s is the unmodulated selection coefficient for deleterious mutations. how the value of s is utilized in the simulation is described above (0 to 1).
* ub is the diploid genome wide rate of beneficial mutations occuring at C-loci (Float).
* hb is the dominance coefficient for beneficial mutations occuring at C-loci (Float).
* sb is the selection coefficient for beneficial mutations occuring at C-loci (0 to 1).
* is is the rate of selfing in the population (0 to 1).
* ir is the rate of recombination, where in each generation there is ir recombination events per individual (Float).
* sto controls how selection fluctuates as described above (1 = bimodal model, 2 = exponential model).
* f is the temporal autocorrelation in selective states over time (0 to 1).
* phi is the probability that selection fluctuates to s\*d0b at F1-loci (0 to 1).
* psi is the probability that selection fluctuates to s\*d1b at F2-loci (0 to 1).
* eps controls the correlation in selective states between F1 and F2 under the bimodal model (Float).
* eta controls the correlation in selective states between F1 and F2 under the exponential model (-1, 0,1).
* d0a controls the fitness effects of deleterious mutation at F1-loci (Float).
* d0b controls the fitness effects of deleterious mutation at F1-loci (Float).
* d1a controls the fitness effects of deleterious mutation at F2-loci (Float).
* d1b controls the fitness effects of deleterious mutation at F2-loci (Float).
* d2c controls the fitness effects of deleterious mutation at C-loci (Float).
* prc is the proportion of loci in the genome under constant selection (0 to 1).
* prf0 is the proporion of F-loci that will be placed in to set F1 (0 to 1).

Below are a few examples of changes that can be made to the sample parameter file above:
1. -prc 0.75 would cause 75% of the genome to contain C-loci and 25% to contain F-loci.  
2. -phi 0.25 would cause the selection coeeificent at F1-loci to be *s d1a* 25% of the time and *s d1b* 75% of the time. 
3. -is 1 would cause the population to be obligately selfing.  
4. -eps 0.25 would cause the correlation in selective states between F1 and F2 to be 1.  
5. -eps -0.25 would cause the correlation in selective states between F1 and F2 to be -1.  
6. -sto 2 would change to the exponential model (values of phi, psi, eps, d0a, d0b, d1a, d1b are all ignored).  

#### OUTPUT
There are three result files outputed from the simulation in the given output directory (all are comma separated files). following the example parameter file above the three output files for the first replicate would be:
1. MySim1_LociGroups.txt  
   This file shows the position of C-loci, F1-loci and F2-loci in the genome.  
 
2. MySim1_EnvFlucData.txt  
   This file gives the frequency and average run length of each environment during the simulation.  
 
3. MySim1.txt  
   This file contains the main results of the simulation. It contains results from generation xmin to mng in intervals of xdat. All relevant parameters are also shown in this output along with information about:  
   * The number of deleterious and beneficial fixations at C-loci and F-loci.
   * The variance in the neutral quantitative trait for calculating effective population size.
   * The mean and variance in fitness among indivduals.
   * The mean and variance number of deleterious and beneficial mutations at C-loci and F-loci across individuals.
