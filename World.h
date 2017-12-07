// This class contains the collection of all populations in the model
// and controls for how many generations the model should run
#include "Popul.h"



// Declaration of the class "theWorld"
class World {
private:

	//vector storing which loci belongs to one of three groups
	//First dimension: fluctuating 1
	//Second dimension: fluctuating 2
	//Third dimension: constant
	vector < vector<int> > vecGroupedLoci;
	
	//Record sequence of environments for DeterministicFluctuation and StochasticFluctuation
	vector <double>	vecEnvType;
	
	//Record sequence of selection coefficients for StochasticFluctuation_Exponential
	vector< vector <double> > vecExpS;
	
	//Store sequence of selection modifiers for all three groups throughout the simulation
	vector< vector<double> > vecDeltaS;
	
	//array to hold current selection modifiers
	double arrDeltaS[3];
	
	//count number of generations that have occured
	int generation;

    //******Eddie MODDED Mar 16 2016
    int fitnessMatches;
    int relativefitnessMatches;
    //******Eddie MODDED Mar 16 2016
    
    //******Eddie MODDED May 4 2016
    //Record data about the selection coefficients occuring during the simulation
    double meanDeltaS[3];
    double varDeltaS[3];
    double meanS[3];
    double varS[3];
    double corS0S1;
	//******Eddie MODDED May 4 2016
	
	//******Eddie MODDED Aug 30 2016
	//3D vector contains all neutral loci values for all individuals at all locus at all allele
	vector< vector< vector<float> > > neutralPop;
	//Same as above, but hold values for the offspring generation
	vector< vector< vector<float> > > neutralOffspring;
	
	//Store mean and var of neutralPop
	double neutralMean[3];
	double neutralVar[3];
	//******Eddie MODDED Aug 30 2016
	
	//*****Eddie MODDED Feb 24 2017
	vector<int> selfLengthPop;
	vector<int> selfLengthOffspring;
	double selfLengthMean;
	double selfLengthVar;
	//*****Eddie MODDED Feb 24 2017
	
	
	double nFocalMut;
	int nFocalLost;
	int nFocalFix;
	int sojournTime;
	int indicateFix;
	double sojournHetero;
        
public:
    
	//Constructor
	World(int simRep, string SaveDir, string SaveFileName);
	
	//Sort loci into three groups (fluctuating 1, fluctuating 2, constant)
	void GenerateLociGroups(vector< vector<int> >& groupedLoci, int simRep, string SaveDir, string SaveFileName);
	
	//Sort loci into two fluctuating groups (Even, Odd)
	void GenerateLociGroups_EvenOdd(vector< vector<int> >& groupedLoci, int simRep, string SaveDir, string SaveFileName);
	
	//Generate sequence of environments under deterministic fluctuation
	void DeterministicFluctuation ( double period, double t0, double t1, double tau, double Gen );
	
	//Generate sequence of environments under stochastic fluctuation between min and max selection
	void StochasticFluctuation ( double f, double phi, double psi, double eps, double Gen, int InitialEnv = -1 );
	
	//Generate sequence of environments under stochastic fluctuation with selection exponentially distributed
	void ConvertEnvToSelection ( double ds0_A, double ds0_B, double ds1_A, double ds1_B, double ds2 );
	
	//Output data about the sequence of environments
	void EnvFluctuationData ( double period, double t0, double t1, double tau, double f, double phi, double psi, double eps, int stochastic, int Round, string SaveDir, string SaveFileName );
	
	//Save relevant data for simulation
	void saveData(Popul& newPop, int whichRounds, string SaveDir, string SaveFileName);
	
	//******Eddie MODDED Mar 16 2016
	void saveSimCheckData(Popul& newPop, int whichRounds, string SaveDir, string SaveFileName);
	void indivFitnessCheck(Popul& newPop, int t);
	void indivRelativeFitnessCheck(Popul& newPop, int t);
	void countMutations(Indiv& newInd, int (&dMut)[3], int (&bMut)[3], int (&dHomo)[3], int (&dHetero)[3], int (&bHomo)[3], int (&bHetero)[3]);
	//******Eddie MODDED Mar 16 2016
	
	//******Eddie MODDED May 3 2016
	//can specify initial s values, default as -1 so will be randomly generated
	void StochasticFluctuation_Exponential ( double f, double beta, double eta, double Gen, double init_s0=-1, double init_s1=-1, double init_s2=-1 );
	//******Eddie MODDED May 3 2016
	
	
	//runModel1: regular simulations with autosomal mutations
	void runModel1(int simRep, string SaveDir, string SaveFileName);
	
	//runModel2: simplified simulation without autosomal mutation, but with random mating and neutrality
	void runModel2(int simRep, string SaveDir, string SaveFileName);
	
	//calc mean and variance of neutral values in neutralPop
	void calcNeutralValues();
	
	//calc mean and variance of selfing length
	void calcSelfLength();
	
	//runModel3_pFix: neutral locus[1] converted into focal locus to track fixation probability
	void runModel3_pFix(int simRep, string SaveDir, string SaveFileName);
	void calcFocalMut();
	void saveData_pFix(Popul& newPop, int simRep, string SaveDir, string SaveFileName);
	
	
};
