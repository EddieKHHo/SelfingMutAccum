#include "World.h"
#include "InputPar.h"



string  MODELFILENAME;

double  NUMOFREPS;
double  NUMOFGENERATIONS;
double  MAXNUMOFGENERATIONS;

int  	L;
int	  	POPULATION;
double  MUTATIONRATE;
double  H;
double  S;
double  K;

double  INITIALSELFING;
double  SELFINGMUTATIONRATE;
double  SIGMAS;
double  FINALSELFINGMUTATIONRATE;
double  FINALSIGMAS;

double  INITIALRECOMBINATION;
double  RECOMBMUTATIONRATE;
double  SIGMAR;

double  STOCHASTIC_FLUCTUATION;

//Deterministic fluctuations
double 	PERIOD;;
double 	T0;
double 	T1;
double 	TAU;

//Stochastic fluctuations
double	F_AC;
double	PHI;
double	PSI;
double	EPSILON;

double  DELTA_S_EVEN_A;
double  DELTA_S_EVEN_B;
double  DELTA_S_ODD_A;
double  DELTA_S_ODD_B;

int		SAVE_DATA_INTERVAL;


//**********Eddie MODDED July 21 2015
double	P_LOCI_GROUP[3];
int		RUNTIME_ERRORS;

double	DELTA_S0_A;
double	DELTA_S0_B;
double	DELTA_S1_A;
double	DELTA_S1_B;
double	DELTA_S2_C;

double	P_CONSTANT;
double	P_FLUC_G0;

double	ENV_FLUC_INDICATOR;
//**********Eddie MODDED July 21 2015

//**********Eddie MODDED Mar 8 2016
double CHECK_FIXATION;
//**********Eddie MODDED Mar 8 2016

double SAVE_DATA_MIN_GEN;
double Q_NEUTRAL_LOCI;

double ETA;

double SELF_LENGTH;

double IGNORE_LOCI;

double PFIX_MODEL;
double PFIX_H;
double PFIX_S;
double PFIX_K;

double U_BEN;
double S_BEN;
double H_BEN;

int main(int argc, char* argv[]) {
    
	//Default parameter values
	MODELFILENAME = "DEFAULT_MODEL_FILE_NAME";
	
    NUMOFREPS = 1;
    NUMOFGENERATIONS = 25000;
    MAXNUMOFGENERATIONS = NUMOFGENERATIONS;

    L = 5000;
    POPULATION = 10000;
    MUTATIONRATE = 1;
    H = 0.25;
    S = 0.05;
    K = 0;
    
    STOCHASTIC_FLUCTUATION = 1;
    
    F_AC = 0.95;
	PHI = 0.5;
	PSI = 0.5;
	EPSILON = 0.25;
	
    PERIOD = 80;
    T0 = 0.5;
    T1 = 0.5;
    TAU = 0;
    
    INITIALSELFING = 1;
    
    SELFINGMUTATIONRATE = 0;
    SIGMAS = 0;
    FINALSELFINGMUTATIONRATE = 0;
    FINALSIGMAS = 0;
    
    INITIALRECOMBINATION = 5;
    RECOMBMUTATIONRATE = 0;
    SIGMAR = 0;

    SAVE_DATA_INTERVAL = 1;    
	
	P_LOCI_GROUP[0] = 0.3;
	P_LOCI_GROUP[1] = 0.3;
	P_LOCI_GROUP[2] = 0.4;
	RUNTIME_ERRORS = 0;
	DELTA_S0_A = 1;
	DELTA_S0_B = 1;
	DELTA_S1_A = 1;
	DELTA_S1_B = 1;
	DELTA_S2_C = 1;
	
	P_CONSTANT = 0;
	P_FLUC_G0 = 0.5;
	RUNTIME_ERRORS = 0;
	
	CHECK_FIXATION = 0;
	
	SAVE_DATA_MIN_GEN = 1;
	Q_NEUTRAL_LOCI = 0;
	ETA = 0;
	
	SELF_LENGTH = 0;
	
	IGNORE_LOCI = 0;
	
	PFIX_MODEL = 0;
	PFIX_H = 0;
	PFIX_S = 0;
	PFIX_K = 0;
	
	U_BEN = 0;
	S_BEN = 0;
	H_BEN = 0;
	//Default parameter values

	
    
    int 						Errors = 0;
    vector<string> 				Keys;
	vector< vector<double> >	Values;
	vector<string> 				ResultFileNames;
	string 						SaveDirectory;
	int 						NumberSimulations = 0;
	int 						NumberSimulationNames = 0;
	
	//Generate ExampleParameterFile.txt  just in case ParameterFileName not specified
	//ExampleParameterFile.txt will NOT allow simulation to start because it does NOT specify a save directory
	WriteExampleParameterFile();
	
	//Default to example parameter file existing in working directory of terminal
	string						ParameterFileName = "ExampleParameterFile.txt";
	
	//Store USER specified parameter file directory/name in ParameterFileName
	if( argc > 1){
		ParameterFileName = argv[1];
	}
	
	//Read parameters, save directory and save names from specified file
	//If file not specified, defaults to ExampleParameterFile.txt in the working directory of terminal    
   	ReadParameterFile (Keys, Values, ResultFileNames, SaveDirectory, NumberSimulations, NumberSimulationNames, ParameterFileName);
    
    //Check initially to make sure Save directory is specified and NumberSimulations == NumberSimulationNames
	//Simulation will NOT continure if any of this fails
	Errors += InitialErrorCheck(SaveDirectory, NumberSimulations, NumberSimulationNames);
    
    
    //Continus if 0 initial errors
	if( Errors == 0 ){
		
		//MakeParameterList takes Values (2D vector) of parameters, creates the parameter lists and store it in ParameterList	
		vector<double> emptyVec;
		vector< vector<double> > ParameterList;
		MakeParameterList(Values,0,emptyVec,ParameterList);
	
		//output ParList on screen
		cout << "Print: ParList" << "\n";
		PrintVecDouble2D(ParameterList);
		cout << "\n\n";
	
		for(int i=0; i < NumberSimulations; i++){	
		
			//Assign MODELFILENAME, which is file name of the results file		
			MODELFILENAME = ResultFileNames[i];
			
			//Assign parameter values in ParameterList to global variables of simulation
			AssignParameterValues (Keys, ParameterList, i, SaveDirectory, MODELFILENAME);
			
			//Check parameters for errors, return number to Errors
			Errors += ParameterErrorCheck(SaveDirectory, MODELFILENAME);
				
			//If one or more errors, simulation will NOT run, check ErrorFile if this occurs
			if(Errors == 0){			
			
				//perform NUMOFREPS replicate simulations
				for(int i = 1; i <= NUMOFREPS; i++) {
					
					if(PFIX_MODEL == 0){
						//Run model in modelWorld
						World modelWorld( i, SaveDirectory, MODELFILENAME );
					
						modelWorld.runModel1( i, SaveDirectory, MODELFILENAME);
						//modelWorld.runModel2( i, SaveDirectory, MODELFILENAME);
					}
					else if(PFIX_MODEL == 1){
						//Run model in modelWorld
						World modelWorld( i, SaveDirectory, MODELFILENAME );
					
						modelWorld.runModel3_pFix( i, SaveDirectory, MODELFILENAME);
					}

				}
				
				cout << "File name of results is: " << MODELFILENAME << "\n";
			}
			else{
				cout << "Simulation stopped for: " << MODELFILENAME << "\n";
				cout << "FAILED parameter error check. " << "\n";
				cout << "Please look at Error File: " << SaveDirectory << MODELFILENAME << "_ParErrorFile.txt" << "\n";
			}
			
			//RESTART Error check for next simulation, regardless if current simulation stopped
			//Don't want parameter errors from one simulation to be inherited by next one, which may pass parameters error check
			Errors = 0;
		}		
	}
	//fail error check, never started simulation
	else{
		cout << "Simulation stopped." << "\n";
		cout << "FAILED initial error check due to unspecified save directory and/or number of simulations not equal to number of simulation names given."<< "\n";
		cout << "Please look at Error File: "<< SaveDirectory << "InitialErrorFile.txt" << "\n";
	}
	
	return 0;
}



