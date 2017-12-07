// This class represents individuals in the model
#include "IncludeFiles.h"
#include "SharedFunctions.h"

//Declare friend class prototype
class Popul;
class World;



//Declaration of the class "indiv"
class Indiv {

private:
	
	//2D vector, first index is locus group 0,1,2; second index is loci containing del/ben mutation within that group
	vector < vector<int> > Cm1;
	vector < vector<int> > Cm2;

	//arrays to store value of modifier loci for recomb and self in Cm1 (index 0) and Cm2 (index 1)
	double modifierRecomb[2];
	double modifierSelfing[2];
	
	//individual fitness, pollen discounted fitness, selfing rate and recombination rate
	long double indivFitness;
	long double indivPDFitness;
	double indivSelfingRate;
	double indivRecombRate;
	
	//arrays to store num del/ben mut in loci groups 0,1,2
	int delMut[3];
	int benMut[3];
	//arrays to store num of del/ben homozygous and heterozygous loci in loci groups 0,1,2
	int delHomo[3];
	int delHetero[3];
	int benHomo[3];
	int benHetero[3];
	
	//ID indicate whether individual was create from selfing(0) or outcrossing(1)
	int ID;
	

public:

	//***Indiv has FOUR constructors for with different purporses
	//   Two are for runModel1, other two are for runModel2	
	//initialize first generation of individuals (not used for rest of simulation)
	Indiv( vector < vector<int> >& groupedLoci, double (&deltaS)[3] );
	
	//sets gamete1,2 recomb rate, self rate, and neutral loci values for a new individual
	Indiv( vector < vector<int> > gamete1, vector < vector<int> > gamete2, double recomb1, double recomb2, double selfing1, double selfing2);
	
	void addInitialMutation( vector < vector<int> >& groupedLoci );
	
	void calcFitness_V3( double (&deltaS)[3] );
	
	void calcSelfingRate();
	
	void calcRecombRate();
	
	//function references to variables that will store info about newly created gamete
	void makeGamete(double& meanRecombPopul, vector< vector<int> >& gamete, double& gameteModifierRecomb, double& gameteModifierSelfing, int (&gameteNeutralAllele)[3]);
	

	//*****Eddie MODDED Aug 1 2016
	Indiv();
	Indiv(double recomb1, double recomb2, double selfing1, double selfing2);
	void makeGamete_New(double& gameteModifierRecomb, double& gameteModifierSelfing, int (&gameteNeutralCm)[3]);
	//*****Eddie MODDED Aug 1 2016
		
	void calcGenotype_V3();
	
	
	//class studyPop and theWorld can access all the information of indiv
	friend class Popul;
	friend class World;
};
