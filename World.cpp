#include "World.h"

//World constructor
//1. Classify loci in to group, write groups into file
//2. Generate temporal sequence of env either deterministically or stochastically
//3. Convert env to selection coefficient
//4. Write env data in to files (has expected and observed values)
//5. Display basic info for what type of environment generated
//6. Display if using neutral loci
//7. Initialize neutralPop is using neutral loci
World::World(int simRep, string SaveDir, string SaveFileName){
	
	if(PFIX_MODEL==1){
		cout << "#--------------------------------------------------#\n";
		cout << "Running World_pFix\n";
		cout << "Will keep track of a focal locus in middle of genome for fix probability\n";
		cout << "Neutral locus[1] is now under selection\n";
		cout << "#--------------------------------------------------#\n";
	}

	//Generate three sets of loci groups    
	GenerateLociGroups(vecGroupedLoci, simRep, SaveDir, SaveFileName);
	
	//***Use ONLY if want G0 to have all Odd mutations and G1 to have all even mutations
	//GenerateLociGroups_EvenOdd(vecGroupedLoci, simRep, SaveDir, SaveFileName);
	
	
	if(STOCHASTIC_FLUCTUATION==0){
		DeterministicFluctuation( PERIOD, T0, T1, TAU, MAXNUMOFGENERATIONS);
	}
	else if (STOCHASTIC_FLUCTUATION==1){
		StochasticFluctuation ( F_AC, PHI, PSI, EPSILON, MAXNUMOFGENERATIONS );
	}
	else if (STOCHASTIC_FLUCTUATION==2){
		double BETA = S/2;
		
		StochasticFluctuation_Exponential(F_AC, BETA, ETA, MAXNUMOFGENERATIONS);
	}
		
	//Convert sequence of environments to sequence of selection coefficients over time		
	ConvertEnvToSelection ( DELTA_S0_A, DELTA_S0_B, DELTA_S1_A, DELTA_S1_B, DELTA_S2_C );
		
	if(STOCHASTIC_FLUCTUATION==0 || STOCHASTIC_FLUCTUATION==1){	
		//Save data of sequence of env changes (allow comparison to expected values)
		EnvFluctuationData ( PERIOD, T0, T1, TAU, F_AC, PHI, PSI, EPSILON, STOCHASTIC_FLUCTUATION, simRep, SaveDir, SaveFileName );
	}
		
	//initial DeltaS values
	arrDeltaS[0] = vecDeltaS[0][0];
	arrDeltaS[1] = vecDeltaS[1][0];
	arrDeltaS[2] = vecDeltaS[2][0];
	
	
	//*****Display Simulation Info
	if(STOCHASTIC_FLUCTUATION == 0){
		cout << "Simulation is DETERMINISTIC.\n";
		cout << meanS[0] << "," << meanS[1] << "," << meanS[2] << "\n";
		cout << varS[0] << "," << varS[1] << "," << varS[2] << "\n";
		cout << "Correlation between s0 and s1: " << calcCorrelation_Double(vecDeltaS[0],vecDeltaS[1]) << "\n";
	}
	else if (STOCHASTIC_FLUCTUATION == 1){
		cout << "Simulation is STOCHASTIC.\n";
		cout << meanS[0] << "," << meanS[1] << "," << meanS[2] << "\n";
		cout << varS[0] << "," << varS[1] << "," << varS[2] << "\n";
		cout << "Correlation between s0 and s1: " << calcCorrelation_Double(vecDeltaS[0],vecDeltaS[1]) << "\n";
	}
	else if (STOCHASTIC_FLUCTUATION == 2){
		cout << "Simulation is STOCHASTIC EXPONENTIAL\n";
		cout << meanS[0] << "," << meanS[1] << "," << meanS[2] << "\n";
		cout << varS[0] << "," << varS[1] << "," << varS[2] << "\n";
		cout << "Correlation between s0 and s1: " << calcCorrelation_Double(vecExpS[0],vecExpS[1]) << "\n";
	}
	
	//determine is tracking quantitative neutral loci
	if(Q_NEUTRAL_LOCI == 0 && PFIX_MODEL==0){
		cout << "Neutral loci NOT IN USE." << "\n";
	}
	else if (Q_NEUTRAL_LOCI == 1 | PFIX_MODEL==1){
		
		// Initialize the random number generator
		const gsl_rng_type * T;
		gsl_rng * r;
		T = gsl_rng_taus;
		r = gsl_rng_alloc (T);
		gsl_rng_default_seed += time(NULL);
		gsl_rng_set (r, gsl_rng_default_seed);
		
		cout << "Neutral loci IN USE." << "\n";
		
		//*****Eddie MODDED Aug, 30, 2016
		//initialize neutralPop with POPULATION,3,2 elements in dimensions 1,2,3 respectively
		neutralPop.resize(POPULATION);
		for(int n = 0; n < POPULATION; n++){
			neutralPop[n].resize(3);
		
			for(int l = 0; l < 3; l++){
				neutralPop[n][l].resize(2);
			
				if(PFIX_MODEL==0){
					//start neutral values as random gaussian variates with a variance of 2*POPULATION (st.dev of sqrt(2*POPULATION))
					neutralPop[n][l][0]=gsl_ran_gaussian(r,sqrt(2*POPULATION));
					neutralPop[n][l][1]=gsl_ran_gaussian(r,sqrt(2*POPULATION));
				}
				else if (PFIX_MODEL==1){
					//start neutral values as 0
					neutralPop[n][l][0]=0;
					neutralPop[n][l][1]=0;
				}
				
			}
		}
		
		//Display initial neutral values
		calcNeutralValues();
		cout << "Initial neutralVar: " << neutralVar[0] << "," << neutralVar[1] << "," << neutralVar[2] << "\n";
		
		//initialize neutralOffspring with POPULATION,3,2 elements in dimensions 1,2,3 respectively
		neutralOffspring.resize(POPULATION);
		for(int n = 0; n < POPULATION; n++){
			neutralOffspring[n].resize(3);
		
			for(int l = 0; l < 3; l++){
				neutralOffspring[n][l].resize(2);
			
				//start all allele neutral value as 0
				neutralOffspring[n][l][0]=0;
				neutralOffspring[n][l][1]=0;
			}
		}
		//*****Eddie MODDED Aug, 30, 2016
		
		//free memory; prevent memory leaks
		gsl_rng_free(r);
	}
	
	//Determine if tracking self length
	if(SELF_LENGTH == 0){
		cout << "NOT tracking consecutive number of generations a lineage undergoes selfing." << "\n";
	}
	else if (SELF_LENGTH == 1){
		cout << "Tracking consecutive number of generations a lineage undergoes selfing." << "\n";
		//initialize vectors that track selfing length
		selfLengthPop.resize(POPULATION);
		selfLengthOffspring.resize(POPULATION);
	}
	
	//Warning about value at IGNORE_LOCI
	if(IGNORE_LOCI == 0){
		cout << "\n#--------------------#\n";
		cout << "IGNORE_LOCI = 0\n";
		cout << "Tracking ALL loci\n";
		cout << "#--------------------#\n";
	}
	else if (IGNORE_LOCI == 1){
		cout << "\n#--------------------#\n";
		cout << "IGNORE_LOCI = 1\n";
		cout << "TRACKING ONLY FLUCTUATING LOCI\n";
		cout << "NOT TRACKING ONLY CONSTANT LOCI\n";
		cout << "#--------------------#\n";
	}
	else if (IGNORE_LOCI == 2){
		cout << "\n" << "#--------------------#\n";
		cout << "IGNORE_LOCI = 2\n";
		cout << "TRACKING ONLY CONSTANT LOCI\n";
		cout << "NOT TRACKING ONLY FLUCTUATING LOCI\n";
		cout << "#--------------------#\n";
	}
}

//Run model
//1. Generate first generation of population (newPop)
//2. Run model for MAXNUMOFGENERATIONS generations
//		Update environment 
//		Calculate fitness and find highest fitness
//		Calc stats and write data
//		Make new generation
//		Remove fixed mutations
//		Update neutral loci values if in use
void World::runModel1(int simRep, string SaveDir, string SaveFileName){
	
	// Initialize the random number generator
	const gsl_rng_type * T;
	gsl_rng * r;
	T = gsl_rng_taus;
	r = gsl_rng_alloc (T);
	gsl_rng_default_seed += time(NULL);
	gsl_rng_set (r, gsl_rng_default_seed);
	
	cout << "Simulation: runModel1" << "\n";
	
	//Create new population (object Popul) containing 2*POPULATION individuals (object Indiv)
	Popul newPop(vecGroupedLoci, arrDeltaS);
	
	//initialize generation as 1
	generation = 1;
	
	cout << "Start simulation: runModel1" << "\n";

	//Run simulation for MAXNUMOFGENERATIONS number of generations
	for(int t = 1; t <= MAXNUMOFGENERATIONS; t++){
	
		//Update selection coefficients
		arrDeltaS[0] = vecDeltaS[0][t];
		arrDeltaS[1] = vecDeltaS[1][t];
		arrDeltaS[2] = vecDeltaS[2][t];
	
		//Calculate fitness of ALL individuals with current selection coefficients
		for(int n = 0; n < POPULATION; n++){
			newPop.thePop[n].calcFitness_V3(arrDeltaS);
		}
		//Find highest fitness and pollen discounted fitness; REQUIRED to get relative fitness
		newPop.findHighestFitness();
		

		//Display data every SAVE_DATA_INTERVAL generations
		//IMPORTANT: SAVE_DATA_MIN_GEN must be >= SAVE_DATA_INTERVAL because header of save file only created at gen SAVE_DATA_INTERVAL
		if(generation%SAVE_DATA_INTERVAL == 0){	
			//Display sim info, no calculations or saving data
			cout << "\n";
			cout << SaveFileName << simRep << "\n";
			cout << "Gen: " << t << "; " << arrDeltaS[0] << "," << arrDeltaS[1] << "," << arrDeltaS[2] << "\n";
			
			//Calc and write data to save file when generation past SAVE_DATA_MIN_GEN AND divisible by SAVE_DATA_INTERVAL
			//All are UNNECESSARY for simulation to run, it is just for saving data
			if(generation >= SAVE_DATA_MIN_GEN){
				//calculate data for saving
				newPop.calcSelfingRate();
				newPop.calcDelMut();
				newPop.calcBenMut();
				newPop.calcPopulFitness();
				
				//display data
				cout << "populFitness: " << newPop.populFitness[0] << "," << newPop.populFitness[1] << "\n";
				cout << "populDelMut[0]: " << newPop.populDel[0][0] << "," << newPop.populDel[0][1] << "," << newPop.populDel[0][2] << " - " << newPop.populTotalDel[0] << "\n";
				cout << "populBenMut[0]: " << newPop.populBen[0][0] << "," << newPop.populBen[0][1] << "," << newPop.populBen[0][2] << " - " << newPop.populTotalBen[0] << "\n";
				cout << "Recomb, Selfing Rate: " << newPop.populRecombRate[0] << "," << newPop.populSelfingRate[0] << "\n";
				cout << "Fixations_Del: " << newPop.fixDel[0] << "," << newPop.fixDel[1] << "," << newPop.fixDel[2] << "\n";
				cout << "Fixations_Ben: " << newPop.fixBen[0] << "," << newPop.fixBen[1] << "," << newPop.fixBen[2] << "\n";
				//cout << "ID_P: " << newPop.ID_P1 << "," << newPop.ID_P2 << "\n";
				
				if(Q_NEUTRAL_LOCI == 1){
					calcNeutralValues();
					cout << "varNeutral: " << neutralVar[0] << "," << neutralVar[1] << "," << neutralVar[2] << "\n";
				}
				if(SELF_LENGTH == 1){
					calcSelfLength();
					cout << "selfLength: " << selfLengthMean << "," << selfLengthVar << "\n";
				}
				
				cout << "\n";
					
				//write data to save file
				saveData(newPop, simRep, SaveDir, SaveFileName);
			
				//write data if running simulation to check fixation algorithm
				if(CHECK_FIXATION > 0){
					saveSimCheckData(newPop, simRep, SaveDir, SaveFileName);
				}
			}
		}
		
		//Make new generation of individuals: selection, recombination, mutation, reproduction
		newPop.makeNewGeneration(vecGroupedLoci);
		//remove fixations
		newPop.removeFixations();
		
		//if using neutral loci, update values
		if(Q_NEUTRAL_LOCI == 1){		
			//*****Update neutral values
			//1) create offspring vector of neutral values based of which parents pass on alleles to next generation
			//2) mutate all neutral values
			//3) neutralPop updated as neutralOffspring
			//4) erase values at neutralOffspring
			for(int n = 0; n < POPULATION; n++){		
				for(int l = 0; l < 3; l++){			
					//neutralOffspring[n][l][] inherit neutralPop neutral values based on indices in neutralGameteIndex[n][l][][]
					neutralOffspring[n][l][0] = neutralPop[ newPop.neutralGameteIndex[n][l][0][0] ][l][ newPop.neutralGameteIndex[n][l][0][1] ];
					neutralOffspring[n][l][1] = neutralPop[ newPop.neutralGameteIndex[n][l][1][0] ][l][ newPop.neutralGameteIndex[n][l][1][1] ];
				
					//mutate neutral values
					neutralOffspring[n][l][0] += gsl_ran_gaussian(r,1.0);
					neutralOffspring[n][l][1] += gsl_ran_gaussian(r,1.0);
				}
			}
		
			//update neutralPop using offspring generation
			neutralPop = neutralOffspring;
			//*****Update neutral values
			
		}
		
		//if tracking self lengths, update selfing lengths
		if(SELF_LENGTH == 1){
			//*****Update selfing lengths
			//newPop.neutralGameteIndex[n][l][0][0] and newPop.neutralGameteIndex[n][l][1][0] indicates which parent offspring n will inherit allele 0 and allele 1
			//	the 'l' doesn't matter for parent identify, its only used for calc neutral loci stuff
			for(int n = 0; n < POPULATION; n++){
				
				int Parent1 = newPop.neutralGameteIndex[n][0][0][0];
				int Parent2 = newPop.neutralGameteIndex[n][0][1][0];
					
				if (Parent1 == Parent2){
					selfLengthMean ++;
					selfLengthOffspring[n] = selfLengthPop[Parent1] + 1;
				}
				else{
					selfLengthOffspring[n] = 0;
				}
			
			}
			
			//update selfLengthPop
			selfLengthPop = selfLengthOffspring;
			
			//reset selfLengthOffspring
			fill(selfLengthOffspring.begin(),selfLengthOffspring.end(),0);
			//*****Update selfing lengths
		}
		
		//increase generation count
		generation++;	
	}

	//free memory; prevent memory leaks
	gsl_rng_free(r);
}



//*****EDDIE MODDED Nov 16 2015
void World::GenerateLociGroups_EvenOdd(vector< vector<int> >& groupedLoci, int simRep, string SaveDir, string SaveFileName){
	
	//initialize groupedLoci with three empty vectors (each one represent one class of loci)
	for(int i = 0; i < 3; i++){
		groupedLoci.push_back( vector<int> () );
	}

	//Loop through loci 1 to L
	for(int i = 1; i <= L; i++){
		//Assign even number loci to group[0]
		if(i%2 == 0){
			groupedLoci[0].push_back(i);
		}
		//Assign odd number loci to group[1]
		else{
			groupedLoci[1].push_back(i);
		}
	}
	
	//sort loci in ascending order within group 0,1,2
	sort(groupedLoci[0].begin(), groupedLoci[0].end());
	sort(groupedLoci[1].begin(), groupedLoci[1].end());

	//display size of each loci group
	cout << groupedLoci[0].size() << "," << groupedLoci[1].size() << "," << groupedLoci[2].size() << "\n";
	
	//*****Record loci in each loci group
	// string to indicate which round we are in (round = replicate simulation)
	stringstream sstr;
	sstr << simRep;
	//file name starts are "model results" and sstr (=whichRounds)
	string fileName(SaveDir);
	fileName.append(SaveFileName);
	fileName.append(sstr.str());
	fileName.append("_LociGroups");
	fileName.append(".txt");

	// open the text file in which the data will be saved
	ofstream modelResults;
	modelResults.open(fileName.c_str(), ios::out);
	
	//loop through all three loci groups
	for(int g = 0; g < 3; g++){
		//print current group
		modelResults		<< "LociGroup" << g << "(" << groupedLoci[g].size() << "):";
	
		//print all loci in group g
		for(int i = 0; i < groupedLoci[g].size(); i++){
			modelResults 	<< groupedLoci[g][i] << ",";
		}
		//next line
		modelResults 		<< "\n";
	}
	//close file	
	modelResults.close();
	//*****Record loci in each loci group
}
//*****EDDIE MODDED Nov 16 2015



//Classify loci in to G0, G1 or G2 (G0 and G1 can fluctuate in selection, G2 always constant selection)
//Proportion P_CONSTANT of the loci belong in G2
//Proportion (1-P_CONSTANT)*P_FLUC_G0 in G0 and (1-P_CONSTANT)*(1-P_FLUC_G0) in G1
//Number of loci in each group is calculated
//Loci are chosen from UNIFORM distribution ~(1,L) (ignoring repeat) and placed in to groups based on MULTINOMIAL distribution
//	If a group is full, it cannot get additional loci even if it is chosen by multinomial distribution
void World::GenerateLociGroups(vector< vector<int> >& groupedLoci, int simRep, string SaveDir, string SaveFileName){
	
	//initialize groupedLoci with three empty vectors (each one represent one class of loci)
	for(int i = 0; i < 3; i++){
		groupedLoci.push_back( vector<int> () );
	}

	//initialize count of loci in each of three groups with 0
	vector <int> countLociInGroup(3,0);
	
	//Variables to create probability bins
	vector <double>		MutProb(3,0);
	vector <double> 	Bin(3,0);
	double 				TempNum = 0;
	
	//Assign probabilities of each locus group occurring
	MutProb[0] = P_LOCI_GROUP[0];
	MutProb[1] = P_LOCI_GROUP[1];
	MutProb[2] = P_LOCI_GROUP[2];
	
	//calc num of loci to assign to each group, total = L
	vector <int> numLociToAssign(3,0);
	numLociToAssign[0] = floor(L * (1 - P_CONSTANT) * P_FLUC_G0);
	numLociToAssign[1] = floor(L * (1 - P_CONSTANT) * (1 - P_FLUC_G0) );
	numLociToAssign[2] = L - numLociToAssign[0] - numLociToAssign[1];
	cout << numLociToAssign[0] << "," << numLociToAssign[1] << "," << numLociToAssign[2] << "\n";
	
	
	//*****Sept 24
	
	//Calc the boundary for bins that will categorize locus group along interval 0 to 1
	for(int i = 0; i < 3; i++){
		TempNum	+= MutProb[i];
		Bin[i] = TempNum;
	}
		
	// Initialize the random number generator
	const gsl_rng_type * T;
	gsl_rng * r;
	T = gsl_rng_taus;
	r = gsl_rng_alloc (T);
	gsl_rng_default_seed += time(NULL);
	gsl_rng_set (r, gsl_rng_default_seed);
	
	//store locus that has been assigned to prevent duplication
	vector<int> vecAssignedLoci;
	//count num of loci assigned to a class
	int numAssignedLoci = 0;
	
	//assign random locus to locus group 0,1,2; continue until all L loci are assigned
	while( numAssignedLoci < L){
		
		//indicate if sample locus is a duplicate within record of sampled loci within vecAssignedLoci
		int duplicateLocus;

		//exit do loop only if tempLocus is not found within vecAssignedLoci (i.e. has never been sampled)
		do{
			//initially assume not duplicated
			duplicateLocus = 0;
			//sample random locus between 1 and L
			int tempLocus = floor(gsl_ran_flat(r, 0.5, L+0.4999) + 0.5);
			
			//check if tempLocus already sampled
			for(int k = 0; k < vecAssignedLoci.size(); k++){
				//duplicateLocus = 1 if tempLocus exist within vecAssignedLoci
				if(tempLocus == vecAssignedLoci[k]){
					duplicateLocus = 1;
				}
			}
			
			//if duplicateLocus == 0 after check, tempLocus is not a duplicate, assign it to a non-full loci group
			if(duplicateLocus == 0){
				//record tempLocus in vecAssignedLoci
				vecAssignedLoci.push_back(tempLocus);
				//increase count of loci assigned
				numAssignedLoci++;
				
				//cout << tempLocus << "," << vecAssignedLoci.size() << "\n";
								
				//indicate whether tempLocus is assigned to a group
				int tempLocusAssigned = 0;
				
				//while loop continues until tempLocusAssigned == 1 (i.e. when it is assigned to a group)
				//	tempLocusAssigned remains 0 ONLY if it is assigned to a certain class that is full
				//		if that occurs, loop will choose a new random uniform variate 
				//Assignment to group depends on probability Bins and require that it has not been filled up yet
				while(tempLocusAssigned == 0){
					//random uniform variate
					double UnifRV = gsl_ran_flat(r, 0.0, 1.0);	
					gsl_rng_default_seed += time(NULL);
					gsl_rng_set (r, gsl_rng_default_seed);
					
					
					//assign to loci group 0,1,2 according to probability bins AND only if they are NOT full
					if( UnifRV <= Bin[0] && countLociInGroup[0] < numLociToAssign[0] ){
						//record tempLocus in mut group 0
						groupedLoci[0].push_back(tempLocus);
						//increase count of loci assigned to group 0
						countLociInGroup[0]++;
						//indicate tempLocus has been assigned; exit while loop
						tempLocusAssigned = 1;
					}
					else if ( UnifRV > Bin[0] && UnifRV <= Bin[1] && countLociInGroup[1] < numLociToAssign[1] ){
						groupedLoci[1].push_back(tempLocus);
						countLociInGroup[1]++;
						tempLocusAssigned = 1;
					}
					else if( countLociInGroup[2] < numLociToAssign[2] ){
						groupedLoci[2].push_back(tempLocus);
						countLociInGroup[2]++;
						tempLocusAssigned = 1;
					}
				}
			}
		}while(duplicateLocus == 1);	
	}
	
	//sort loci in ascending order within group 0,1,2
	sort(groupedLoci[0].begin(), groupedLoci[0].end());
	sort(groupedLoci[1].begin(), groupedLoci[1].end());
	sort(groupedLoci[2].begin(), groupedLoci[2].end());

	cout << groupedLoci[0].size() << "," << groupedLoci[1].size() << "," << groupedLoci[2].size() << "\n";
	
	
	//*****Record loci in each loci group
	// string to indicate which round we are in (round = replicate simulation)
	stringstream sstr;
	sstr << simRep;
	//file name starts are "model results" and sstr (=whichRounds)
	string fileName(SaveDir);
	fileName.append(SaveFileName);
	fileName.append(sstr.str());
	fileName.append("_LociGroups");
	fileName.append(".txt");

	// open the text file in which the data will be saved
	ofstream modelResults;
	modelResults.open(fileName.c_str(), ios::out);
	
	//loop through all three loci groups
	for(int g = 0; g < 3; g++){
		//print current group
		modelResults		<< "LociGroup" << g << "(" << groupedLoci[g].size() << "):";
	
		//print all loci in group g
		for(int i = 0; i < groupedLoci[g].size(); i++){
			modelResults 	<< groupedLoci[g][i] << ",";
		}
		//next line
		modelResults 		<< "\n";
	}
	//close file	
	modelResults.close();
	//*****Record loci in each loci group
	
	//free memory; prevent memory leaks
	gsl_rng_free(r);
}



//Generate sequence of environments based on deterministic algorithm (PERIOD, T0, T1, TAU)
//Assume there is only 4 types of environments {Even/Odd} = {A/A, A/B, B/A, B/B}
void World::DeterministicFluctuation ( double period, double t0, double t1, double tau, double Gen ){
	
	//clear vector
	vecEnvType.clear();
	
	//initialize vecEnvType with 0
	for( int i = 0; i <= Gen; i++){
		vecEnvType.push_back(0);
	}
	
	//i0, i1 = 1 if in env A and = 0 if in env B
	double i0;
	double i1;
	
	double genCount=1;
	
	for(int i = 0; i <= Gen; i++){
        //Even in env A from genCount = 1 to PERIOD*T0
        if(genCount <= period*t0){
            i0 = 1;
        }
        else{
            i0 = 0;
        }
        //Odd in env A from genCount = (PERIOD*T0*TAU) to ((PERIOD*T0*TAU)+PERIOD*T1)
        if(genCount > (period*t0*tau) && genCount <= ((period*t0*tau)+period*t1)){
            i1 = 1;
        }
        else{
            i1 = 0;
        }
        
        //Record type of env depending on status of i0 and i1
        if( i0 == 1 && i1 == 1){
        	vecEnvType[i] = 1;
        }
        else if ( i0 == 1 && i1 == 0){
        	vecEnvType[i] = 2;
        }
        else if ( i0 == 0 && i1 == 1){
        	vecEnvType[i] = 3;
        }
        else{
        	vecEnvType[i] = 4;
        }
        
        //when genCount==PERIOD, reset it to 1 (AND skip ++)to start cycle again
        if(genCount == period){
            genCount = 1;
        }
        else{
            genCount++;
        }
    }
    
}

//Generate sequence of environments based on stochastic algorithm (f, phi and psi)
//Assume there is only 4 types of environments {Even/Odd} = {A/A, A/B, B/A, B/B}
//Assume vector referenced to vecEnvType is initialized with Gen number of elements using: vector <double> SomeVector(Gen, 0)
//	i.e. size of vecEnvType should be same as Gen
//Optional to parameterize InitialEnv; default InitialEnv = -1 (written in theWorld.h) and will cause random initial env
void World::StochasticFluctuation ( double f, double phi, double psi, double eps, double Gen, int InitialEnv ){

	//NumEnv needs to be a const, so I can initialize arrays with size [NumEnv];
	const int			NumEnv = 4;
	vector <double>		EnvProb(NumEnv,0);
	vector <double> 	Bin(NumEnv,0);
	double 				TempNum = 0;

	//Assign probabilities of each environment occurring
	EnvProb[0] = (phi)*(psi) + eps;
	EnvProb[1] = (phi)*(1-psi) - eps;
	EnvProb[2] = (1-phi)*(psi) - eps;
	EnvProb[3] = (1-phi)*(1-psi) + eps;
	
	//Calculate the boundary for the bins that will categorize environments along interval 0 to 1
	for(int i = 0; i < NumEnv; i++){
		TempNum	+= EnvProb[i];
		Bin[i] = TempNum;
	}
	
	// Initialize the random number generator
	const gsl_rng_type * T;
	gsl_rng * r;
	T = gsl_rng_taus;
	r = gsl_rng_alloc (T);
	gsl_rng_default_seed += time(NULL);
	gsl_rng_set (r, gsl_rng_default_seed);

	
	//Initialize CurrentEnv
	double 	CurrentEnv;
	
	//if InitialEnv between 1 and 4, CurrentEnv start as InitialEnv
	if(InitialEnv >= 1 &&  InitialEnv <= 4){
		CurrentEnv = InitialEnv;
		
		cout << "Using given initial env" << "\n";
	}
	//otherwise, random generate CurrentEnv
	else{
		double 	UnifRV = gsl_ran_flat(r, 0.0, 1.0);	
		gsl_rng_default_seed += time(NULL);
		gsl_rng_set (r, gsl_rng_default_seed);

		if( UnifRV <= Bin[0] ){
			CurrentEnv = 1;
		}
		else if ( UnifRV > Bin[0] && UnifRV <= Bin[1] ){
			CurrentEnv = 2;
		}
		else if ( UnifRV > Bin[1] && UnifRV <= Bin[2] ){
			CurrentEnv = 3;
		}
		else{
			CurrentEnv = 4;
		}
		cout << "Random generate initial env" << "\n";
	}
	
	//clear vector
	vecEnvType.clear();
	
	//initialize vecEnvType with 0
	for( int i = 0; i <= Gen; i++){
		vecEnvType.push_back(0);
	}
	
	//Generate sequence of environments for Gen cycles from index 0 to Gen
	//	vecEnvType[0] = CurrentEnv (initialized above) ONLY if it doesn't change after one generation
	//	This is done so that the "user-defined initial CurrentEnv" is actually the LAST env of a previous sequence subject to change in the new sequence being generated
	//	If we use vecEnvType[0] = "user-defined initial CurrentEnv = LAST env of a previous sequence", there will be no possible chance of fluctuation from LAST env of a previous sequence to FIRST env of new sequence (which is incorrect)
	//Simulate multinomial sampling by picking uniform RV (between 0 and 1)
	//Bins correspond to size of each environment, env determined by which bin uniform RV falls into
	for(int i = 0; i <= Gen; i++){

		//Fluctuate = {0,1} means {remain the same, fluctuate}
		//Using 1-f because f denotes amount of autocorrelation/no change; we want prob of success for changing env = 1-f
		double Fluctuate = gsl_ran_binomial (r, 1-f, 1);
		gsl_rng_default_seed += time(NULL);
		gsl_rng_set (r, gsl_rng_default_seed);
		
		//Env fluctuate when Fluctuate == 1
		//Update CurrentEnv to result of fluctuation (possible to fluctuate to same env)
		if( Fluctuate == 1){
		
			double UnifRV = gsl_ran_flat(r, 0.0, 1.0);	
	
			if( UnifRV <= Bin[0] ){
				CurrentEnv = 1;
			}
			else if ( UnifRV > Bin[0] && UnifRV <= Bin[1] ){
				CurrentEnv = 2;
			}
			else if ( UnifRV > Bin[1] && UnifRV <= Bin[2] ){
				CurrentEnv = 3;
			}
			else{
				CurrentEnv= 4;
			}
		}
		
		//cout << UnifRV << "," << EnvType[i] << "\n";
		
		//Save EnvType
		vecEnvType[i] = CurrentEnv;
	}
	
	//free memory; prevent memory leaks
	gsl_rng_free(r);			
}



//Generate sequence of environments based on stochastic algorithm with selection being exponentially distributed
void World::StochasticFluctuation_Exponential ( double f, double beta, double eta, double Gen, double init_s0, double init_s1, double init_s2 ){

	// Initialize the random number generator
	const gsl_rng_type * T;
	gsl_rng * r;
	T = gsl_rng_taus;
	r = gsl_rng_alloc (T);
	gsl_rng_default_seed += time(NULL);
	gsl_rng_set (r, gsl_rng_default_seed);

	//temp_s will store s values for current generation
	double temp_s[3];
	
	//if InitialEnv between 0 and 1, CurrentEnv start as InitialEnv
	if(init_s0 >= 0 &&  init_s0 <= 1 && init_s1 >= 0 &&  init_s1 <= 1 && init_s2 >= 0 &&  init_s2 <= 1){
		temp_s[0] = init_s0;
		temp_s[1] = init_s1;
		temp_s[2] = init_s2;
		
		cout << "Using given initial s: " << temp_s[0] << "," << temp_s[1] << "," << temp_s[2] << "\n";
	}
	//otherwise, random generate CurrentEnv; exponential dist with mean beta
	else{
		
		//do loop to ensure all temp_s values less than 1 (prob of >1 is extremely small)
		do{
			if(eta== 0){
				//both s0 and s1 random sample from exp dist
				temp_s[0]= gsl_ran_exponential (r,beta);
				temp_s[1]= gsl_ran_exponential (r,beta);
				temp_s[2]= beta;
			}
			else if(eta== 1){
				temp_s[0]= gsl_ran_exponential (r,beta);
				//s1 exactly same as s0
				temp_s[1]= temp_s[0];
				temp_s[2]= beta;
			}
			else if(eta== -1){
				temp_s[0]= gsl_ran_exponential (r,beta);
				//make s1 take on value at 1-quantile of s0
				double quatile_s0 = gsl_cdf_exponential_Q(temp_s[0], beta);
				temp_s[1] = gsl_cdf_exponential_Qinv(1-quatile_s0, beta);
				temp_s[2]= beta;
			}
			else{
				cout << "UNRECOGNIZED eta" << "\n";
			}
		}while(temp_s[0]>1 || temp_s[1]>1 || temp_s[2]>1);
				
		cout << "Random generate initial s" << "\n";
	}
	
	//clear vecExpS
	vecExpS.clear();
	//initialize vecExpS with 0
	for(int g = 0; g < 3; g++){
		vecExpS.push_back( vector<double> ());
	}
	for( int i = 0; i <= Gen; i++){
		vecExpS[0].push_back(0);
		vecExpS[1].push_back(0);
		vecExpS[2].push_back(0);
	}

	
	//Generate sequence of environments for Gen cycles from index 0 to Gen
	//Fluctuation occur with frequency 1-f; then choose CurrentEnv based on linear probability
	for(int i = 0; i <= Gen; i++){

		//Fluctuate = {0,1} means {remain the same, fluctuate}
		//Using 1-f because f denotes amount of autocorrelation/no change; we want prob of success for changing env = 1-f
		double Fluctuate = gsl_ran_binomial (r, 1-f, 1);
		gsl_rng_default_seed += time(NULL);
		gsl_rng_set (r, gsl_rng_default_seed);
		
		//Env fluctuate when Fluctuate == 1
		//Update CurrentEnv to result of fluctuation (possible to fluctuate to same env)
		if( Fluctuate == 1){
		
			//do loop to ensure all temp_s values less than 1 (prob of >1 is extremely small)
			do{
				if(eta== 0){
					//both s0 and s1 random sample from exp dist
					temp_s[0]= gsl_ran_exponential (r,beta);
					temp_s[1]= gsl_ran_exponential (r,beta);
					temp_s[2]= beta;
				}
				else if(eta== 1){
					temp_s[0]= gsl_ran_exponential (r,beta);
					//s1 exactly same as s0
					temp_s[1]= temp_s[0];
					temp_s[2]= beta;
				}
				else if(eta== -1){
					temp_s[0]= gsl_ran_exponential (r,beta);
					//make s1 take on value at 1-quantile of s0
					double quatile_s0 = gsl_cdf_exponential_Q(temp_s[0], beta);
					temp_s[1] = gsl_cdf_exponential_Qinv(1-quatile_s0, beta);
					temp_s[2]= beta;
				}
				else{
					cout << "UNRECOGNIZED eta" << "\n";
				}
			}while(temp_s[0]>1 || temp_s[1]>1 || temp_s[2]>1);
		}
		
				
		//Save vecExpS
		for(int g = 0; g < 3; g++){
			vecExpS[g][i]=temp_s[g];
		}
	}
	
	//free memory; prevent memory leaks
	gsl_rng_free (r);	

}



//Input selection modifiers for Even and Odd loci for every generation based on the type of environment
//Assumed vectors referenced to DS_Even DS_Odd initialized with Gen elements
void World::ConvertEnvToSelection ( double ds0_A, double ds0_B, double ds1_A, double ds1_B, double ds2 ){
	
	double Gen;

	if(STOCHASTIC_FLUCTUATION==0 || STOCHASTIC_FLUCTUATION==1){	
		Gen = vecEnvType.size();
	}
	else if(STOCHASTIC_FLUCTUATION==2){
		Gen = vecExpS[0].size();
	}
	
	
	//local variable to record selection coefficients (S*deltaS)
	vector< vector<double> > vecS;
	
	//clear vector
	vecDeltaS.clear();
	
	//initialize vecDeltaS with 0 for all mut groups
	for( int g = 0; g < 3; g++){
		vecDeltaS.push_back( vector<double> () );
		vecS.push_back(vector<double> () );
	}
	for( int i = 0; i < Gen; i++){
		vecDeltaS[0].push_back(0);
		vecDeltaS[1].push_back(0);
		vecDeltaS[2].push_back(0);
		
		vecS[0].push_back(0);
		vecS[1].push_back(0);
		vecS[2].push_back(0);
	}
	if(STOCHASTIC_FLUCTUATION==0 || STOCHASTIC_FLUCTUATION==1){
		//Since Gen = vecEnvType.size(), we only go from i to Gen-1 rather from i to Gen
		for(int i = 0; i < Gen; i++){
		
			//group 0 and 1 mutations may have variable selective strength over time
			if(vecEnvType[i] == 1){
				vecDeltaS[0][i] = ds0_A;
				vecDeltaS[1][i] = ds1_A;
			}
			else if (vecEnvType[i] == 2){
				vecDeltaS[0][i] = ds0_A;
				vecDeltaS[1][i] = ds1_B;
			}
			else if (vecEnvType[i] == 3){
				vecDeltaS[0][i] = ds0_B;
				vecDeltaS[1][i] = ds1_A;
			}
			else if (vecEnvType[i] == 4){
				vecDeltaS[0][i] = ds0_B;
				vecDeltaS[1][i] = ds1_B;
			}
			else{
				cout << "UNKNOWN ENVIRONMENT TYPE !" << "\n";
			}
		
			//groups 2 mutations are under constant selection
			vecDeltaS[2][i] = ds2;
			
			//save selection values
			vecS[0][i] = vecDeltaS[0][i]*S;
			vecS[1][i] = vecDeltaS[1][i]*S;
			vecS[2][i] = vecDeltaS[2][i]*S;
		}
	}
	else if(STOCHASTIC_FLUCTUATION==2){
			
		//Since Gen = vecEnvType.size(), we only go from i to Gen-1 rather from i to Gen
		for(int i = 0; i < Gen; i++){			
			//vecDeltaS = vecExpS/S so that when calculating fitness, s = S*vecDeltaS = S*vecExpS/S = vecExpS
			vecDeltaS[0][i] = vecExpS[0][i]/S;
			vecDeltaS[1][i] = vecExpS[1][i]/S;
			vecDeltaS[2][i] = vecExpS[2][i]/S;
			
			//save selection values
			vecS[0][i] = vecDeltaS[0][i]*S;
			vecS[1][i] = vecDeltaS[1][i]*S;
			vecS[2][i] = vecDeltaS[2][i]*S;
		}
	}
	
	
	//calculate mean vecDeltaS and mean vecS
	for( int g = 0; g < 3; g++){
		calcMeanVar_Double(vecDeltaS[g],meanDeltaS[g],varDeltaS[g]);
		calcMeanVar_Double(vecS[g],meanS[g],varS[g]);
	}
	//calculate correlation between s0 and s1
	corS0S1=calcCorrelation_Double(vecS[0],vecS[1]);
}


//Record env fluctuation data from vecEnvType (frequency and run length of different environments)
//Check vecEnvType results match expectations based on parameters (deterministic or stochastic)
//ASSUMES only 4 types of env {A/A, A/B, B/A, B/B}
void World::EnvFluctuationData ( double period, double even_A, double odd_A, double tau, double f, double phi, double psi, double eps, int stochastic, int Round, string SaveDir, string SaveFileName ){
	
	//NumEnv needs to be a const, so I can initialize arrays with size [NumEnv];
	const int			NumEnv = 4;
	
	//EnvCount will record number of gen of a certain env occurring
	//EnvRunLength will record run length of a certain env; EnvRuns count number of runs for a certain env
	//EnvCurrent and RunLength are temporary variables for counting
	vector <double> 	EnvCount(NumEnv,0);
	vector <double>		EnvRunLength(NumEnv,0);
	vector <double>		EnvRuns(NumEnv,0);
	double				EnvCurrent;
	double				RunLength;
	
	//*****EDDIE MODDED Nov 16 2015
	
	//Record mean DeltaS per env type (there should be no variance in selection modifier within a particular env type)
	//First dimension is env type, second dimension is mut type
	double meanDeltaS[4][3] = {0};
	
	//*****EDDIE MODDED Nov 16 2015
	
	//Determine num gen in vecEnvType
	double Gen = vecEnvType.size();
	//Since Gen = vecEnvType.size(), we only go from i to Gen-1 rather from i to Gen
	for(int i = 0; i < Gen; i++){
		
		//Initialize EnvCurrent, RunLength and EnvRuns
		if( i == 0 ){
			EnvCurrent = vecEnvType[i];
			RunLength = 0;
		}
		
	
		//Record frequency of each type of env; regardless if it has changed or not
		if( vecEnvType[i] == 1 ){
			EnvCount[0]++;
			meanDeltaS[0][0]+= vecDeltaS[0][i];
			meanDeltaS[0][1]+= vecDeltaS[1][i];
			meanDeltaS[0][2]+= vecDeltaS[2][i];
		}
		else if ( vecEnvType[i] == 2 ){
			EnvCount[1]++;
			meanDeltaS[1][0]+= vecDeltaS[0][i];
			meanDeltaS[1][1]+= vecDeltaS[1][i];
			meanDeltaS[1][2]+= vecDeltaS[2][i];
		}
		else if ( vecEnvType[i] == 3 ){
			EnvCount[2]++;
			meanDeltaS[2][0]+= vecDeltaS[0][i];
			meanDeltaS[2][1]+= vecDeltaS[1][i];
			meanDeltaS[2][2]+= vecDeltaS[2][i];
		}
		else {
			EnvCount[3]++;
			meanDeltaS[3][0]+= vecDeltaS[0][i];
			meanDeltaS[3][1]+= vecDeltaS[1][i];
			meanDeltaS[3][2]+= vecDeltaS[2][i];
		}
		
		//If env has changed since previous OR last Gen (i.e. env hasn't changed, but sim ending)
		//Record RunLength of environment, EnvRuns increase by 1; Set EnvCurrent to new env and reset RunLength
		//Slight inaccuracies occur because we may need record a run length of an env BEFORE it ends due to limited # gen
		if( EnvCurrent != vecEnvType[i] || i == (Gen-1) ){
			
			//Add 1 to RunLength to count the last gen when at the end of simulation
			if( i == (Gen-1) ){
				RunLength ++;
			}
		
			if( EnvCurrent == 1 ){
				EnvRunLength[0] += RunLength;
				EnvRuns[0] ++;
			}
			else if ( EnvCurrent == 2 ){
				EnvRunLength[1] += RunLength;
				EnvRuns[1] ++;
			}
			else if ( EnvCurrent == 3 ){
				EnvRunLength[2] += RunLength;
				EnvRuns[2] ++;
			}
			else {
				EnvRunLength[3] += RunLength;
				EnvRuns[3] ++;
			}
			//update current env, reset run length counter	
			EnvCurrent = vecEnvType[i];
			RunLength = 0;
		}
		//RunLength + 1
		RunLength ++;
	}
	
	//This section calculates the observed freq and run length of different env AND the expected values given the various parameters
	//It will print out expected and observed value in to EnvFlucData.txt; these should match very closely, otherwise, something wrong with simulation
	vector <double> ExpEnvFreq(NumEnv,0);
	vector <double> ObsEnvFreq(NumEnv,0);
	vector <double> ExpEnvRunLength(NumEnv,0);
	vector <double> ObsEnvRunLength(NumEnv,0);
	
	//Calc OBSERVED frequency of each env
	ObsEnvFreq[0] = EnvCount[0]/Gen;
	ObsEnvFreq[1] = EnvCount[1]/Gen;
	ObsEnvFreq[2] = EnvCount[2]/Gen;
	ObsEnvFreq[3] = EnvCount[3]/Gen;
	
	//Calc OBSERVED run length of each env
	ObsEnvRunLength[0] = EnvRunLength[0]/EnvRuns[0];
	ObsEnvRunLength[1] = EnvRunLength[1]/EnvRuns[1];
	ObsEnvRunLength[2] = EnvRunLength[2]/EnvRuns[2];
	ObsEnvRunLength[3] = EnvRunLength[3]/EnvRuns[3];
	
	//Calc OBSERVED mean DeltaS for each env
	meanDeltaS[0][0]/=EnvCount[0];
	meanDeltaS[0][1]/=EnvCount[0];
	meanDeltaS[0][2]/=EnvCount[0];
	meanDeltaS[1][0]/=EnvCount[1];
	meanDeltaS[1][1]/=EnvCount[1];
	meanDeltaS[1][2]/=EnvCount[1];
	meanDeltaS[2][0]/=EnvCount[2];
	meanDeltaS[2][1]/=EnvCount[2];
	meanDeltaS[2][2]/=EnvCount[2];
	meanDeltaS[3][0]/=EnvCount[3];
	meanDeltaS[3][1]/=EnvCount[3];
	meanDeltaS[3][2]/=EnvCount[3];
	
	// string to indicate which round we are in (round = replicate simulation)
	stringstream sstr;
	sstr << Round;
	//file name starts as "FileName" and sstr (=whichRounds)
	string OutputFileName(SaveDir);
	OutputFileName.append(SaveFileName);
	OutputFileName.append(sstr.str());
	OutputFileName.append("_EnvFlucData");
	OutputFileName.append(".txt");

	//ios::out simply sets the file as an output file
	//add .c_str() to make string FILENAME because .open only takes const char* as parameter
	ofstream OutputFile;
	OutputFile.open (OutputFileName.c_str(), ios::out);
	
	
	//For deterministic fluctuations
	if(stochastic == 0){
		//Case for tau <= 1 and tau > 1 is qualitatively different, need different calculations
		if( tau >= 0 && tau <= 1){
			//Calc expected frequency of each env
			ExpEnvFreq[0] = even_A * (1 - tau);
			ExpEnvFreq[1] = even_A * tau;
			ExpEnvFreq[2] = odd_A - even_A * (1 - tau);
			ExpEnvFreq[3] = 1 - even_A * tau - odd_A;
		
			//Calc expected run length of each env
			ExpEnvRunLength[0] = ExpEnvFreq[0] * period;
			ExpEnvRunLength[1] = ExpEnvFreq[1] * period;
			ExpEnvRunLength[2] = ExpEnvFreq[2] * period;
			ExpEnvRunLength[3] = ExpEnvFreq[3] * period;
		}
		//when tau > 1, there are two phases with a B/B env, these are the two terms of ExpEnvFreq[3]
		//hence, ExpEnvRunLength[3] needs to use te average of those two phases to get exp run length of B/B env
		else if ( tau > 1 ){
			//Calc expected frequency of each env
			ExpEnvFreq[0] = 0;
			ExpEnvFreq[1] = even_A;
			ExpEnvFreq[2] = odd_A;
			ExpEnvFreq[3] = ( even_A * (tau - 1) ) + (1 - even_A * tau - odd_A);
		
			//Calc expected run length of each env
			ExpEnvRunLength[0] = ExpEnvFreq[0] * period;
			ExpEnvRunLength[1] = ExpEnvFreq[1] * period;
			ExpEnvRunLength[2] = ExpEnvFreq[2] * period;
			ExpEnvRunLength[3] = (( even_A * (tau - 1) ) + (1 - even_A * tau - odd_A))/2 * period;
		
		
		}
		else{	
			OutputFile << "ERROR: tau < 0" << "\n";
		}	
		
	}
	//For stochastic fluctuations
	else if (stochastic == 1){
		
		//Calc expected frequency of each env
		ExpEnvFreq[0] = (phi)*(psi) + eps;
		ExpEnvFreq[1] = (phi)*(1-psi) - eps;
		ExpEnvFreq[2] = (1-phi)*(psi) - eps;
		ExpEnvFreq[3] = (1-phi)*(1-psi) + eps;
		
		//Calc expected run length of each env
		ExpEnvRunLength[0] = 1/((1-f)*(1-ExpEnvFreq[0]));
		ExpEnvRunLength[1] = 1/((1-f)*(1-ExpEnvFreq[1]));
		ExpEnvRunLength[2] = 1/((1-f)*(1-ExpEnvFreq[2]));
		ExpEnvRunLength[3] = 1/((1-f)*(1-ExpEnvFreq[3]));
		
	}
	else{
		OutputFile << "Env fluctuation algorithm unspecified (stochastic != 0 or 1)" << "\n";
	}
	
	//Write env data to file
		OutputFile << "ModelFileName" << "," << "Total_Gen" << "," << "Stochastic" << "," << "Env" << "," << "ExpEnvFreq" << "," << "ObsEnvFreq" << "," << "ExpEnvRunLength" << "," << "ObsEnvRunLength" << "," << "meanDeltaS[0]" << "," << "meanDeltaS[1]" << "," << "meanDeltaS[2]" <<"\n";
	for(int i = 0; i < NumEnv; i++){
		OutputFile << SaveFileName << "," << Gen << "," << stochastic << "," << (i+1) << "," << ExpEnvFreq[i] << "," << ObsEnvFreq[i] << "," << ExpEnvRunLength[i] << "," << ObsEnvRunLength[i] << "," << meanDeltaS[i][0] << "," << meanDeltaS[i][1] << "," << meanDeltaS[i][2] << "\n";
	}
	
	OutputFile.close();

}


//function for writing data to text file
void World::saveData(Popul& newPop, int simRep, string SaveDir, string SaveFileName){

	//Saving info for simulation
	// string to indicate which round we are in (round = replicate simulation)
	stringstream sstr;
	sstr << simRep;
	//file name starts are "model results" and sstr (=whichRounds)
	string fileName(SaveDir);
	fileName.append(SaveFileName);
	fileName.append(sstr.str());
	fileName.append(".txt");
	
	//First line for modelResults is made ONLY generation == SAVE_DATA_INTERVAL (i.e. the first time we save data)
	if (generation == SAVE_DATA_MIN_GEN){
		// open the text file in which the data will be saved
		ofstream modelResults;
		modelResults.open(fileName.c_str(), ios::out);
		
		
		//**********Eddie MODDED; Change list of variables in ouput
		modelResults	<< "ModelFileName,";
		modelResults	<< "Q_NEUTRAL_LOCI,SELF_LENGTH,IGNORE_LOCI,";
		modelResults	<< "SAVE_DATA_INTERVAL,SAVE_DATA_MIN_GEN,";
		modelResults	<< "Generation,";
		modelResults	<< "Loci,k,N,U,h,s,";
		modelResults	<< "Ub,hb,sb,";
		modelResults	<< "DELTA_S0_A,DELTA_S0_B,DELTA_S1_A,DELTA_S1_B,DELTA_S2_C,";
		modelResults	<< "P_CONSTANT,P_FLUC_G0,";
		modelResults	<< "ENV_FLUC_INDICATOR,";
		modelResults    << "PERIOD,T0,T1,TAU,";
		modelResults	<< "Stochastic_Fluctuation,F_AC,PHI,PSI,EPSILON,ETA,";
		modelResults	<< "Init_Recomb,U_Recomb,Size_Recomb,";
		modelResults    << "Init_Self,Init_USelf,Init_SizeSelf,Final_USelf,Final_SizeSelf,";
		modelResults	<< "arrDeltaS[0],arrDeltaS[1],arrDeltaS[2],";

		modelResults    << "MeanRecomb,VarRecomb,";
		modelResults    << "MeanSelf,VarSelf,";
		modelResults    << "MeanFitness,VarFitness,";
		modelResults    << "MeanDel0,MeanDel1,MeanDel2,";
		modelResults    << "VarDel0,VarDel1,VarDel2,";
		modelResults	<< "FixDel0,FixDel1,FixDel2,";
		modelResults	<< "VarNeutral0,VarNeutral1,VarNeutral2,";
		modelResults	<< "meanDeltaS0,meanDeltaS1,meanDeltaS2,";
		modelResults	<< "varDeltaS0,varDeltaS1,varDeltaS2,";
		modelResults	<< "meanS0,meanS1,meanS2,";
		modelResults	<< "varS0,varS1,varS2,";
		modelResults	<< "corS0S1,";
		modelResults	<< "selfLengthMean,selfLengthVar,";
		modelResults	<< "FixBen0,FixBen1,FixBen2,";		
		modelResults    << "MeanBen0,MeanBen1,MeanBen2,";
		modelResults    << "VarBen0,VarBen1,VarBen2,";
		modelResults    << "\n";
		
		modelResults.close();
	}
	
	// open the text file in which the data will be saved
	ofstream modelResults;
	modelResults.open(fileName.c_str(), ios::out | ios::app);
	
	//print data is modelResults can be opened
	if(modelResults.is_open()){
		// Write info modelResults
		modelResults	<< SaveFileName << ",";
		modelResults	<< Q_NEUTRAL_LOCI << "," << SELF_LENGTH << "," << IGNORE_LOCI << ",";
		modelResults	<< SAVE_DATA_INTERVAL << "," << SAVE_DATA_MIN_GEN << ",";
		modelResults	<< generation << ",";
		modelResults	<< L << "," << K << "," << POPULATION << "," << MUTATIONRATE << "," << H << "," << S << ",";
		modelResults	<< U_BEN << "," << H_BEN << "," << S_BEN << ",";
		modelResults	<< DELTA_S0_A << "," << DELTA_S0_B << "," << DELTA_S1_A << "," << DELTA_S1_B << "," << DELTA_S2_C << ",";
		modelResults	<< P_CONSTANT << "," << P_FLUC_G0 << ",";
		modelResults	<< ENV_FLUC_INDICATOR << ",";
		modelResults    << PERIOD << "," << T0 << "," << T1 << "," << TAU << ",";
		modelResults	<< STOCHASTIC_FLUCTUATION << "," << F_AC << "," << PHI << "," << PSI << "," << EPSILON << "," << ETA << ",";
		modelResults	<< INITIALRECOMBINATION << "," << RECOMBMUTATIONRATE << "," << SIGMAR << ",";
		modelResults    << INITIALSELFING << "," << SELFINGMUTATIONRATE << "," <<  SIGMAS << "," <<  FINALSELFINGMUTATIONRATE << "," << FINALSIGMAS << ",";
		modelResults	<< arrDeltaS[0] << "," << arrDeltaS[1] << "," << arrDeltaS[2] << ",";
		
		modelResults    << newPop.populRecombRate[0] << "," << newPop.populRecombRate[1] << ",";
		modelResults    << newPop.populSelfingRate[0] << "," << newPop.populSelfingRate[1] << ",";
		modelResults	<< newPop.populFitness[0] << "," << newPop.populFitness[1] << ",";
		modelResults	<< newPop.populDel[0][0] << "," << newPop.populDel[0][1] << "," << newPop.populDel[0][2] << ",";
		modelResults	<< newPop.populDel[1][0] << "," << newPop.populDel[1][1] << "," << newPop.populDel[1][2] << ",";
		modelResults    << newPop.fixDel[0] << "," << newPop.fixDel[1] << "," << newPop.fixDel[2] << ",";
		modelResults    << neutralVar[0] << "," << neutralVar[1] << "," << neutralVar[2] << ",";
		modelResults	<< meanDeltaS[0] << "," << meanDeltaS[1] << "," << meanDeltaS[2] << ",";
		modelResults	<< varDeltaS[0] << "," << varDeltaS[1] << "," << varDeltaS[2] << ",";
		modelResults	<< meanS[0] << "," << meanS[1] << "," << meanS[2] << ",";
		modelResults	<< varS[0] << "," << varS[1] << "," << varS[2] << ",";
		modelResults	<< corS0S1 << ",";
		modelResults	<< selfLengthMean << "," << selfLengthVar << ",";
		modelResults	<< newPop.fixBen[0] << "," << newPop.fixBen[1] << "," << newPop.fixBen[2] << ",";
		modelResults	<< newPop.populBen[0][0] << "," << newPop.populBen[0][1] << "," << newPop.populBen[0][2] << ",";
		modelResults	<< newPop.populBen[1][0] << "," << newPop.populBen[1][1] << "," << newPop.populBen[1][2] << ",";
		modelResults    << "\n";	
	}
	else{
		cout << "CANNOT OPEN modelResults FILE!!!!!" << "\n";
	}
	modelResults.close();
}


//function for writing data to text file
void World::saveSimCheckData(Popul& newPop, int simRep, string SaveDir, string SaveFileName){

	//Saving info for simulation
	// string to indicate which round we are in (round = replicate simulation)
	stringstream sstr;
	sstr << simRep;
	//file name starts are "model results" and sstr (=whichRounds)
	string fileName(SaveDir);
	fileName.append(SaveFileName);
	fileName.append(sstr.str());
	fileName.append("_SimCheck");
	fileName.append(".txt");
	
	//First line for modelResults is made ONLY generation == SAVE_DATA_INTERVAL (i.e. the first time we save data)
	if (generation == SAVE_DATA_MIN_GEN){
		// open the text file in which the data will be saved
		ofstream modelResults;
		modelResults.open(fileName.c_str(), ios::out);
		
		
		//**********Eddie MODDED; Change list of variables in ouput
		modelResults	<< "ModelFileName,";
		modelResults	<< "Generation,";
		modelResults	<< "Loci,k,N,U,h,s,";
		modelResults	<< "DELTA_S0_A,DELTA_S0_B,DELTA_S1_A,DELTA_S1_B,DELTA_S2_C,";
		modelResults	<< "P_CONSTANT,P_FLUC_G0,";
		modelResults	<< "ENV_FLUC_INDICATOR,";
		modelResults    << "PERIOD,T0,T1,TAU,";
		modelResults	<< "Stochastic_Fluctuation,F_AC,PHI,PSI,EPSILON,";
		modelResults	<< "Init_Recomb,U_Recomb,Size_Recomb,";
		modelResults    << "Init_Self,Init_USelf,Init_SizeSelf,Final_USelf,Final_SizeSelf,";
		modelResults	<< "arrDeltaS[0],arrDeltaS[1],arrDeltaS[2],";

		modelResults    << "vecFixationSize,fixationMatches,";
		modelResults	<< "FixDel0,FixDel1,FixDel2,";
		modelResults    << "\n";
		
		modelResults.close();
	}
	
	// open the text file in which the data will be saved
	ofstream modelResults;
	modelResults.open(fileName.c_str(), ios::out | ios::app);
	
	//print data is modelResults can be opened
	if(modelResults.is_open()){
		// Write info modelResults
		modelResults	<< SaveFileName << ",";
		modelResults	<< generation << ",";
		modelResults	<< L << "," << K << "," << POPULATION << "," << MUTATIONRATE << "," << H << "," << S << ",";
		modelResults	<< DELTA_S0_A << "," << DELTA_S0_B << "," << DELTA_S1_A << "," << DELTA_S1_B << "," << DELTA_S2_C << ",";
		modelResults	<< P_CONSTANT << "," << P_FLUC_G0 << ",";
		modelResults	<< ENV_FLUC_INDICATOR << ",";
		modelResults    << PERIOD << "," << T0 << "," << T1 << "," << TAU << ",";
		modelResults	<< STOCHASTIC_FLUCTUATION << "," << F_AC << "," << PHI << "," << PSI << "," << EPSILON << ",";
		modelResults	<< INITIALRECOMBINATION << "," << RECOMBMUTATIONRATE << "," << SIGMAR << ",";
		modelResults    << INITIALSELFING << "," << SELFINGMUTATIONRATE << "," <<  SIGMAS << "," <<  FINALSELFINGMUTATIONRATE << "," << FINALSIGMAS << ",";
		modelResults	<< arrDeltaS[0] << "," << arrDeltaS[1] << "," << arrDeltaS[2] << ",";
		
		modelResults    << newPop.vecFixationSize << "," << newPop.fixationMatches << ",";
		modelResults    << newPop.fixDel[0] << "," << newPop.fixDel[1] << "," << newPop.fixDel[2] << ",";
		
		modelResults    << "\n";	
	}
	else{
		cout << "CANNOT OPEN modelResults FILE!!!!!" << "\n";
	}
	modelResults.close();
}



//****Eddie MODDED Mar 16 2016
//Recalculate ABSOLUTE fitness of individuals and check it is same as population
void World::indivFitnessCheck(Popul& newPop, int t){
	//initialize fitnessMatches
	fitnessMatches = 0;
	
	//variable to hold number of mutations and homo, hetero loci
	int dMut[3], bMut[3];
	int dHomo[3], dHetero[3];
	int bHomo[3], bHetero[3];

	//loop through all individuals
	for(int n = 0; n < POPULATION; n++){
		//Count number of mut and homo hetero loci
		countMutations(newPop.thePop[n], dMut, bMut, dHomo, dHetero, bHomo, bHetero);
		
		double checkFitness;
		checkFitness = 	1 * pow( (1 - S * vecDeltaS[0][t]), dHomo[0] ) * pow( (1 - S * H * vecDeltaS[0][t]), dHetero[0] ) *
						pow( (1 - S * vecDeltaS[1][t]), dHomo[1] ) * pow( (1 - S * H * vecDeltaS[1][t]), dHetero[1] ) *
						pow( (1 - S * vecDeltaS[2][t]), dHomo[2] ) * pow( (1 - S * H * vecDeltaS[2][t]), dHetero[2] );
		
		if( abs(newPop.thePop[n].indivFitness - checkFitness) < 0.0001){
			fitnessMatches++;
		}
		else{
			//cout << "Ind " << n << ": " << newPop.thePop[n].indivFitness << "!=" << checkFitness << "\n";
		}
	}

}

//Recalculate RELATIVE fitness of individuals and check it is same as population
void World::indivRelativeFitnessCheck(Popul& newPop, int t){

	//variable to hold number of mutations and homo, hetero loci
	int dMut[3], bMut[3];
	int dHomo[3], dHetero[3];
	int bHomo[3], bHetero[3];
	
	vector< double > checkFitness;
	double maxFitness;
	
	//loop through all individuals
	for(int n = 0; n < POPULATION; n++){
		//Count number of mut and homo hetero loci
		countMutations(newPop.thePop[n], dMut, bMut, dHomo, dHetero, bHomo, bHetero);
		
		checkFitness.push_back(	1 * pow( (1 - S * vecDeltaS[0][t]), dHomo[0] ) * pow( (1 - S * H * vecDeltaS[0][t]), dHetero[0] ) *
									pow( (1 - S * vecDeltaS[1][t]), dHomo[1] ) * pow( (1 - S * H * vecDeltaS[1][t]), dHetero[1] ) *
									pow( (1 - S * vecDeltaS[2][t]), dHomo[2] ) * pow( (1 - S * H * vecDeltaS[2][t]), dHetero[2] ) );
		if(n == 0){
			maxFitness = checkFitness[n];
		}
		else{
			if(maxFitness < checkFitness[n]){
				maxFitness = checkFitness[n];
			}			
		}
	}

	//initialize fitnessMatches
	relativefitnessMatches = 0;	
	//Check if relative fitness values match
	for(int n = 0; n < POPULATION; n++){
	
		if( abs(newPop.thePop[n].indivFitness/newPop.populHighestFitness - checkFitness[n]/maxFitness) < 0.0001){
			relativefitnessMatches++;
		}
		else{
			//cout << "Ind " << n << ": " << newPop.thePop[n].indivFitness/newPop.populHighestFitness << "!=" << checkFitness[n]/maxFitness << "\n";
		}
	}
}

//Count number homozygous and heterozygous loci for each loci group given an individual (class Indiv)
//Method is different and less efficient than calcFitness_V3()
//input as reference: Indiv newInd, arrays to contain counts of del and ben homozygous and heterozygous for each loci group (array has 3 elements, one for each loci group)
void World::countMutations(Indiv& newInd, int (&dMut)[3], int (&bMut)[3], int (&dHomo)[3], int (&dHetero)[3], int (&bHomo)[3], int (&bHetero)[3]){
	
	//loop through all i loci groups
	for(int i = 0; i < 3; i++){
	
		//*****Count number of del and ben loci in Cm1 + Cm2 
		//count number of Del and Ben loci
		dMut[i] = 0;
		bMut[i] = 0;
	
		//Only count ben mut if PROB_BEN >0
		if(PROB_BEN > 0){
			for( vector<int>::iterator itI = newInd.Cm1[i].begin(); itI < newInd.Cm1[i].end(); itI++ ){
				if( *itI < 0 ){
					bMut[i]++;
				}
			}
			for( vector<int>::iterator itI = newInd.Cm2[i].begin(); itI < newInd.Cm2[i].end(); itI++ ){
				if( *itI < 0 ){
					bMut[i]++;
				}
			}
		}
		
		//num del loci is total subtract num ben loci
		dMut[i] = newInd.Cm1[i].size() + newInd.Cm2[i].size() - bMut[i];
		//*****Count number of del and ben loci in Cm1 + Cm2 
		
		dHomo[i] = 0;
		bHomo[i] = 0;
		dHetero[i] = 0;
		bHetero[i] = 0;
	
		//*****Find number of homozygous and heterozygous loci (del and ben)
		//Use different and less efficient method than calcFitness_V3; it goes through all of Cm1 and Cm2 to find matches
		//loop through all j mutations Cm1
		for(int j = 0; j < newInd.Cm1[i].size(); j++){
			//set focal mutation in Cm1
			int mutCm1=newInd.Cm1[i][j];
			
			//loop through all k mutations in Cm2 until find match or reach end
			for(int k = 0; k < newInd.Cm2[i].size(); k++){
				//set mutation in Cm2 to compare to Cm1
				int mutCm2=newInd.Cm2[i][k];
				
				//increase homozygous count if match
				if(abs(mutCm1) == abs(mutCm2)){
					if(mutCm1 > 0 && mutCm2 > 0){	
						dHomo[i]++;
					}
					else if( mutCm1 < 0 && mutCm2 < 0 ){	
						bHomo[i]++;
					}
					
					//break out of k loop if found match
					break;
				}
			}
		}
		
		//calculate num of heterozygous loci after looking through all mutations
		dHetero[i] = dMut[i] - 2 * dHomo[i];
		bHetero[i] = bMut[i] - 2 * bHomo[i];
		//*****Find number of homozygous and heterozygous loci (del and ben)
	}
}
//****Eddie MODDED Mar 16 2016





//*****Eddie MODDED Aug 2016
//Model without tracking autosome mutations and with no selection; random mating and neutral
void World::runModel2(int simRep, string SaveDir, string SaveFileName){
	
	//ensures Q_NEUTRAL_LOCI = 1 since runModel2 only tracking neutral loci
	Q_NEUTRAL_LOCI = 1;
	
	// Initialize the random number generator
	const gsl_rng_type * T;
	gsl_rng * r;
	T = gsl_rng_taus;
	r = gsl_rng_alloc (T);
	gsl_rng_default_seed += time(NULL);
	gsl_rng_set (r, gsl_rng_default_seed);
	
	cout << "Simulation: runModel2" << "\n";
	
	//*****Create new population (object Popul) containing 2*POPULATION individuals (object Indiv)
	Popul newPop;
	
	//initialize generation as 1
	generation = 1;
	
	cout << "Start simulation: runModel2" << "\n";

	//Run simulation for MAXNUMOFGENERATIONS number of generations
	for(int t = 1; t <= MAXNUMOFGENERATIONS; t++){
	
		//Update selection coefficients
		arrDeltaS[0] = vecDeltaS[0][t];
		arrDeltaS[1] = vecDeltaS[1][t];
		arrDeltaS[2] = vecDeltaS[2][t];
	
		
		//Display data every SAVE_DATA_INTERVAL generations
		//IMPORTANT: SAVE_DATA_MIN_GEN must be >= SAVE_DATA_INTERVAL because header of save file only created at gen SAVE_DATA_INTERVAL
		if(generation%SAVE_DATA_INTERVAL == 0){	
			//Display sim info, no calculations or saving data
			cout << "\n";
			cout << SaveFileName << simRep << "\n";
			cout << "Gen: " << t << "; " << arrDeltaS[0] << "," << arrDeltaS[1] << "," << arrDeltaS[2] << "\n";
			
			//Calc and write data to save file when generation past SAVE_DATA_MIN_GEN AND divisible by SAVE_DATA_INTERVAL
			//All are UNNECESSARY for simulation to run, it is just for saving data
			if(generation >= SAVE_DATA_MIN_GEN){
				//calculate data for saving
				newPop.calcSelfingRate();
				newPop.calcDelMut();
				newPop.calcPopulFitness();
				
				calcNeutralValues();
				
				//display data
				cout << "Recomb, Selfing Rate: " << newPop.populRecombRate[0] << "," << newPop.populSelfingRate[0] << "\n";
				cout << "ID_P: " << newPop.ID_P1 << "," << newPop.ID_P2 << "\n";
				
				cout << "varNeutral: " << neutralVar[0] << "," << neutralVar[1] << "," << neutralVar[2] << "\n";
				cout << "\n";
					
				//write data to save file
				saveData(newPop, simRep, SaveDir, SaveFileName);
			
			}
		}
		
		//Make new generation of individuals: selection, recombination, mutation, reproduction
		newPop.makeNewGeneration_New();

		
		//*****Update neutral values
		//1) create offspring vector of neutral values based of which parents pass on alleles to next generation
		//2) mutate all neutral values
		//3) neutralPop updated as neutralOffspring
		//4) erase values at neutralOffspring
		for(int n = 0; n < POPULATION; n++){		
			for(int l = 0; l < 3; l++){
				//neutralOffspring[n][l][] inherit neutralPop neutral values based on indices in neutralGameteIndex[n][l][][]
				neutralOffspring[n][l][0] = neutralPop[ newPop.neutralGameteIndex[n][l][0][0] ][l][ newPop.neutralGameteIndex[n][l][0][1] ];
				neutralOffspring[n][l][1] = neutralPop[ newPop.neutralGameteIndex[n][l][1][0] ][l][ newPop.neutralGameteIndex[n][l][1][1] ];
				
				//mutate neutral values
				neutralOffspring[n][l][0] += gsl_ran_gaussian(r,1.0);
				neutralOffspring[n][l][1] += gsl_ran_gaussian(r,1.0);
			}
		}
		
		//update neutralPop using offspring generation
		neutralPop = neutralOffspring;
		//*****Update neutral values
		
		//increase generation count
		generation++;	
	}
	
	//free memory; prevent memory leaks
	gsl_rng_free(r);
}







//Calc mean and var num of neutral loci value
void World::calcNeutralValues(){
	//initialize with zeros
	for(int l = 0; l < 3; l++){
		neutralMean[l] = 0;
		neutralVar[l] = 0;
	}
	
	//store sum squares
	double neutralSS[3] = {};
		
	//Get values from all inds
	for(int n = 0; n < POPULATION; n++){
		
		//loop through all three neutral loci
		for(int l = 0; l < 3; l++){
			//store value in neutralLoci[j][0] and neutralLoci[j][1] of ind i
			double temp0 = neutralPop[n][l][0];
			double temp1 = neutralPop[n][l][1];
			
			neutralMean[l] += temp0;
			neutralMean[l] += temp1;
			
			neutralSS[l] += pow(temp0, 2.0);
			neutralSS[l] += pow(temp1, 2.0);

		}
	}
		
	//Calc mean and var for each neutral loci with respect to number of loci (2N)
	for(int l = 0; l < 3; l++){
		neutralMean[l] /= 2*POPULATION;
		neutralVar[l] = neutralSS[l]/(2*POPULATION) - pow( neutralMean[l], 2.0 );
	}

}


void World::calcSelfLength(){
	//initialize mean and var and sum of squares
	selfLengthMean = 0;
	selfLengthVar = 0;
	double selfLengthSS = 0;
	
	//Get values from all inds
	for(int n = 0; n < POPULATION; n++){
		selfLengthMean += selfLengthPop[n];
		selfLengthSS += pow(selfLengthPop[n], 2.0);
	}
	
	//Calc mean and var
	selfLengthMean /= POPULATION;
	selfLengthVar = selfLengthSS/POPULATION - pow(selfLengthMean, 2.0);

}










//Run model
//1. Generate first generation of population (newPop)
//2. Run model for MAXNUMOFGENERATIONS generations
//		Update environment 
//		Calculate fitness and find highest fitness
//		Calc stats and write data
//		Make new generation
//		Remove fixed mutations
//		Update neutral loci values if in use
//Using neutral locus[1] as the focal selected locus (no longer neutral)
//	neutral locus [0] and [2] still neutral
void World::runModel3_pFix(int simRep, string SaveDir, string SaveFileName){
	
	
		
	// Initialize the random number generator
	const gsl_rng_type * T;
	gsl_rng * r;
	T = gsl_rng_taus;
	r = gsl_rng_alloc (T);
	gsl_rng_default_seed += time(NULL);
	gsl_rng_set (r, gsl_rng_default_seed);	

	//Initialize variables
	int timerAddMut = 0;
	int addNewMut = 0;
	sojournTime = 0;
	sojournHetero = 0;
	nFocalMut = 0;
	nFocalLost = 0;
	nFocalFix = 0;

	cout << "Simulation: runModel3_pFix" << "\n";
	
	//Create new population (object Popul) containing 2*POPULATION individuals (object Indiv)
	Popul newPop(vecGroupedLoci, arrDeltaS);
	
	//initialize generation as 1
	generation = 1;
	
	cout << "Start simulation: runModel3_pFix" << "\n";
	//Run simulation for MAXNUMOFGENERATIONS number of generations
	for(int t = 1; t <= MAXNUMOFGENERATIONS; ){
		
		
		//-----Add a SINGLE mut to focal locus (neutral locus[1]) if mut extinct-----
		//Mutation added after offspring created (so technically appears after the mutation step, but its ok
		if(nFocalMut == 0 && timerAddMut <= 0 && generation >= 5000){
			
			//indicate new mutation added
			addNewMut = 1;
			
			//pick random individual (1 to POPULATION)
			int indexInd = floor(gsl_ran_flat(r, -0.5, (POPULATION-1)+0.4999) + 0.5);
			//pick random Cm (0 or 1)
			int indexCm = floor(gsl_ran_flat(r, -0.5, 1.4999) + 0.5);
						
			//initialize neutral locus[1] with mutation
			neutralPop[indexInd][1][indexCm] = 1;
		}

		
		//-----Update selection coefficients-----
		arrDeltaS[0] = vecDeltaS[0][t];
		arrDeltaS[1] = vecDeltaS[1][t];
		arrDeltaS[2] = vecDeltaS[2][t];
	
		//-----Calculate fitness of ALL individuals with current selection coefficients-----
		for(int n = 0; n < POPULATION; n++){
			newPop.thePop[n].calcFitness_V3(arrDeltaS);
		}
		
		//---Modify fitness based on mut at focal locus
		for(int n = 0; n < POPULATION; n++){
			//determine how many mut at focal locus for ind n
			int nFocalMutAtInd = neutralPop[n][1][0]+neutralPop[n][1][1];
			
			//heterozygous at focal locus
			if(nFocalMutAtInd == 1){
				newPop.thePop[n].indivFitness *= (1 - PFIX_S * PFIX_H);
			}
			//homozygous mutant at focal locus
			else if(nFocalMutAtInd == 2){
				newPop.thePop[n].indivFitness *= (1 - PFIX_S);
			}
		}
		
		//---Find highest fitness and pollen discounted fitness; REQUIRED to get relative fitness
		newPop.findHighestFitness();
		

		//Display data every SAVE_DATA_INTERVAL generations
		//IMPORTANT: SAVE_DATA_MIN_GEN must be >= SAVE_DATA_INTERVAL because header of save file only created at gen SAVE_DATA_INTERVAL
		if(generation%SAVE_DATA_INTERVAL == 0){	
			//Display sim info, no calculations or saving data
			cout << "\n";
			cout << SaveFileName << simRep << "\n";
			cout << "Gen: " << t << "; " << arrDeltaS[0] << "," << arrDeltaS[1] << "," << arrDeltaS[2] << "\n";
			
			//Calc and write data to save file when generation past SAVE_DATA_MIN_GEN AND divisible by SAVE_DATA_INTERVAL
			//All are UNNECESSARY for simulation to run, it is just for saving data
			if(generation >= SAVE_DATA_MIN_GEN){
				//calculate data for saving
				newPop.calcSelfingRate();
				newPop.calcDelMut();
				newPop.calcPopulFitness();
				
				//display data
				cout << "populFitness: " << newPop.populFitness[0] << "," << newPop.populFitness[1] << "\n";
				cout << "populDelMut[0]: " << newPop.populDel[0][0] << "," << newPop.populDel[0][1] << "," << newPop.populDel[0][2] << " - " << newPop.populTotalDel[0] << "\n";
				cout << "Recomb, Selfing Rate: " << newPop.populRecombRate[0] << "," << newPop.populSelfingRate[0] << "\n";
				cout << "Fixations: " << newPop.fixDel[0] << "," << newPop.fixDel[1] << "," << newPop.fixDel[2] << "\n";
				
				//calc neutral values (locus[1] is focal locus, so not quantitative anymore)
				calcNeutralValues();
				cout << "varNeutral: " << neutralVar[0] << "," << neutralVar[1] << "," << neutralVar[2] << "\n";
				cout << "\n";
					
				//write data to save file
				saveData(newPop, simRep, SaveDir, SaveFileName);
			
			}
		}
		
		//-----Make new generation of individuals: selection, recombination, mutation, reproduction-----
		newPop.makeNewGeneration(vecGroupedLoci);
		//-----remove fixations-----
		newPop.removeFixations();
		

		//-----Update neutral values-----
		//1) create offspring vector of neutral values based of which parents pass on alleles to next generation
		//2) mutate all neutral values
		//3) neutralPop updated as neutralOffspring
		//neutral locus[1] is selected on, but inheritence can be treated the same as neutral alleles (except will not mutate like the neutral loci)
		for(int n = 0; n < POPULATION; n++){		
			for(int l = 0; l < 3; l++){			
				//neutralOffspring[n][l][] inherit neutralPop neutral values based on indices in neutralGameteIndex[n][l][][]
				neutralOffspring[n][l][0] = neutralPop[ newPop.neutralGameteIndex[n][l][0][0] ][l][ newPop.neutralGameteIndex[n][l][0][1] ];
				neutralOffspring[n][l][1] = neutralPop[ newPop.neutralGameteIndex[n][l][1][0] ][l][ newPop.neutralGameteIndex[n][l][1][1] ];
			}
			
			//mutate neutral values for locus 0 and 2 (does not alter locus 1, which is focal locus for pFix)
			neutralOffspring[n][0][0] += gsl_ran_gaussian(r,1.0);
			neutralOffspring[n][0][1] += gsl_ran_gaussian(r,1.0);
			neutralOffspring[n][2][0] += gsl_ran_gaussian(r,1.0);
			neutralOffspring[n][2][1] += gsl_ran_gaussian(r,1.0);
		}
	
		//update neutralPop using offspring generation
		neutralPop = neutralOffspring;
	
		//increase generation count
		generation++;
		t++;
		
		//decrease timer count
		timerAddMut--;
		
		//-----Count how many mut at focal locus
		calcFocalMut();
		
		//-----If new mut added, then track it
		if(addNewMut==1){
			
			sojournTime++;
			sojournHetero += 2*(nFocalMut/(2*POPULATION))*(1-nFocalMut/(2*POPULATION));
		
			//if mutation fixed, record and reset
			if(nFocalMut == 2*POPULATION){
			
				//record fixation
				nFocalFix += 1;
				indicateFix = 1;
				cout << "Mutation fixed\n";
				cout << "Fix, lost, sojourn: " << nFocalFix << "," << nFocalLost << "," << sojournTime <<"\n";
				saveData_pFix(newPop, simRep, SaveDir, SaveFileName);
			
				//reset focal locus
				for(int n = 0; n < POPULATION; n++){
					neutralPop[n][1][0] = 0;
					neutralPop[n][1][1] = 0;
				}
				
				//reset variables
				timerAddMut = 0;
				addNewMut = 0;
				sojournTime = 0;
				sojournHetero = 0;
			}
			//else if mutation lost
			else if(nFocalMut == 0){
			
				//record lost
				nFocalLost +=1;
				indicateFix = 0;
				cout << "Mutation lost\n";
				cout << "Fix, lost, sojourn: " << nFocalFix << "," << nFocalLost << "," << sojournTime <<"\n";
				saveData_pFix(newPop, simRep, SaveDir, SaveFileName);
				
				//reset variables
				timerAddMut = 0;
				addNewMut = 0;
				sojournTime = 0;
				sojournHetero = 0;
		
			}
		}
		
		//break out of for loop if obtain enough fixations
		if( nFocalFix == PFIX_K){
			cout << "DONE after this many mutations: " << nFocalFix + nFocalLost << "\n";
			cout << "Fix, lost: " << nFocalFix << "," << nFocalLost << "\n";
			break;
			
		}
		
		//----If num of max generations almost past, regenerate sequence of environments
		if( (MAXNUMOFGENERATIONS-t) < 0 ){	
			
			cout << "#--------------------------------------------------#\n";
			cout << "Regenerating environmental sequence\n";	
			cout << "#--------------------------------------------------#\n";

			if (STOCHASTIC_FLUCTUATION==1){
				double initEnv = vecEnvType[t];
				StochasticFluctuation ( F_AC, PHI, PSI, EPSILON, MAXNUMOFGENERATIONS, initEnv );
			}
			else if (STOCHASTIC_FLUCTUATION==2){
				double BETA = S/2;
				double inits0 = vecDeltaS[0][t];
				double inits1 = vecDeltaS[0][t];
				double inits2 = vecDeltaS[0][t];
				StochasticFluctuation_Exponential(F_AC, BETA, ETA, MAXNUMOFGENERATIONS, inits0, inits1, inits2);
			}
		
			//Convert sequence of environments to sequence of selection coefficients over time		
			ConvertEnvToSelection ( DELTA_S0_A, DELTA_S0_B, DELTA_S1_A, DELTA_S1_B, DELTA_S2_C );
			
			//reset time
			t = 1;		
		}
	}

	//free memory; prevent memory leaks
	gsl_rng_free(r);
	
}


//update number of mutations at focal locus
void World::calcFocalMut(){

	nFocalMut = 0;
	
	//count number of mut at focal locus (two alleles per individual)
	for(int n = 0; n < POPULATION; n++){
		nFocalMut += neutralPop[n][1][0];
		nFocalMut += neutralPop[n][1][1];
	}

}



//function for writing data to text file
void World::saveData_pFix(Popul& newPop, int simRep, string SaveDir, string SaveFileName){

	//Saving info for simulation
	// string to indicate which round we are in (round = replicate simulation)
	stringstream sstr;
	sstr << simRep;
	//file name starts are "model results" and sstr (=whichRounds)
	string fileName(SaveDir);
	fileName.append(SaveFileName);
	fileName.append(sstr.str());
	fileName.append(".pFix.txt");
	
	//First line for modelResults is made ONLY for first mutational loss or fixation
	if ( (nFocalLost+nFocalFix) == 1){
		// open the text file in which the data will be saved
		ofstream modelResults;
		modelResults.open(fileName.c_str(), ios::out);
		
		
		//**********Eddie MODDED; Change list of variables in ouput
		modelResults	<< "ModelFileName,";
		modelResults	<< "Q_NEUTRAL_LOCI,SELF_LENGTH,IGNORE_LOCI,";
		modelResults	<< "SAVE_DATA_INTERVAL,SAVE_DATA_MIN_GEN,";
		modelResults	<< "Generation,";
		modelResults	<< "Loci,k,N,U,h,s,";
		modelResults	<< "DELTA_S0_A,DELTA_S0_B,DELTA_S1_A,DELTA_S1_B,DELTA_S2_C,";
		modelResults	<< "P_CONSTANT,P_FLUC_G0,";
		modelResults	<< "ENV_FLUC_INDICATOR,";
		modelResults    << "PERIOD,T0,T1,TAU,";
		modelResults	<< "Stochastic_Fluctuation,F_AC,PHI,PSI,EPSILON,ETA,";
		modelResults	<< "Init_Recomb,U_Recomb,Size_Recomb,";
		modelResults    << "Init_Self,Init_USelf,Init_SizeSelf,Final_USelf,Final_SizeSelf,";
		modelResults	<< "arrDeltaS[0],arrDeltaS[1],arrDeltaS[2],";

		modelResults    << "indicateFix,nFocalFix,nFocalLost,";
		modelResults    << "sojournTime,sojournHetero";
		modelResults    << "\n";
		
		modelResults.close();
	}
	
	// open the text file in which the data will be saved
	ofstream modelResults;
	modelResults.open(fileName.c_str(), ios::out | ios::app);
	
	//print data is modelResults can be opened
	if(modelResults.is_open()){
		// Write info modelResults
		modelResults	<< SaveFileName << ",";
		modelResults	<< Q_NEUTRAL_LOCI << "," << SELF_LENGTH << "," << IGNORE_LOCI << ",";
		modelResults	<< SAVE_DATA_INTERVAL << "," << SAVE_DATA_MIN_GEN << ",";
		modelResults	<< generation << ",";
		modelResults	<< L << "," << K << "," << POPULATION << "," << MUTATIONRATE << "," << H << "," << S << ",";
		modelResults	<< DELTA_S0_A << "," << DELTA_S0_B << "," << DELTA_S1_A << "," << DELTA_S1_B << "," << DELTA_S2_C << ",";
		modelResults	<< P_CONSTANT << "," << P_FLUC_G0 << ",";
		modelResults	<< ENV_FLUC_INDICATOR << ",";
		modelResults    << PERIOD << "," << T0 << "," << T1 << "," << TAU << ",";
		modelResults	<< STOCHASTIC_FLUCTUATION << "," << F_AC << "," << PHI << "," << PSI << "," << EPSILON << "," << ETA << ",";
		modelResults	<< INITIALRECOMBINATION << "," << RECOMBMUTATIONRATE << "," << SIGMAR << ",";
		modelResults    << INITIALSELFING << "," << SELFINGMUTATIONRATE << "," <<  SIGMAS << "," <<  FINALSELFINGMUTATIONRATE << "," << FINALSIGMAS << ",";
		modelResults	<< arrDeltaS[0] << "," << arrDeltaS[1] << "," << arrDeltaS[2] << ",";
		
		modelResults	<< indicateFix << "," << nFocalFix << "," << nFocalLost << ",";
		modelResults	<< sojournTime << "," << sojournHetero;
		modelResults    << "\n";	
	}
	else{
		cout << "CANNOT OPEN modelResults FILE!!!!!" << "\n";
	}
	modelResults.close();
}






/*
void World::GenerateLociGroups_MiddleClump(vector< vector<int> >& groupedLoci, int simRep, string SaveDir, string SaveFileName){
	
	//initialize groupedLoci with three empty vectors (each one represent one class of loci)
	for(int i = 0; i < 3; i++){
		groupedLoci.push_back( vector<int> () );
	}
	
	//calc num of loci to assign to each group, total = L
	vector <int> numLociToAssign(3,0);
	numLociToAssign[0] = floor(L * (1 - P_CONSTANT) * P_FLUC_G0);
	numLociToAssign[1] = floor(L * (1 - P_CONSTANT) * (1 - P_FLUC_G0) );
	numLociToAssign[2] = L - numLociToAssign[0] - numLociToAssign[1];
	cout << numLociToAssign[0] << "," << numLociToAssign[1] << "," << numLociToAssign[2] << "\n";
	
	
	numFlank = 
	numCenter = 

	//Loop through loci 1 to L
	for(int i = 1; i <= L; i++){
		//Assign even number loci to group[0]
		if(i%2 == 0){
			groupedLoci[0].push_back(i);
		}
		//Assign odd number loci to group[1]
		else{
			groupedLoci[1].push_back(i);
		}
	}
	
	//sort loci in ascending order within group 0,1,2
	sort(groupedLoci[0].begin(), groupedLoci[0].end());
	sort(groupedLoci[1].begin(), groupedLoci[1].end());

	//display size of each loci group
	cout << groupedLoci[0].size() << "," << groupedLoci[1].size() << "," << groupedLoci[2].size() << "\n";
	
	//*****Record loci in each loci group
	// string to indicate which round we are in (round = replicate simulation)
	stringstream sstr;
	sstr << simRep;
	//file name starts are "model results" and sstr (=whichRounds)
	string fileName(SaveDir);
	fileName.append(SaveFileName);
	fileName.append(sstr.str());
	fileName.append("_LociGroups");
	fileName.append(".txt");

	// open the text file in which the data will be saved
	ofstream modelResults;
	modelResults.open(fileName.c_str(), ios::out);
	
	//loop through all three loci groups
	for(int g = 0; g < 3; g++){
		//print current group
		modelResults		<< "LociGroup" << g << "(" << groupedLoci[g].size() << "):";
	
		//print all loci in group g
		for(int i = 0; i < groupedLoci[g].size(); i++){
			modelResults 	<< groupedLoci[g][i] << ",";
		}
		//next line
		modelResults 		<< "\n";
	}
	//close file	
	modelResults.close();
	//*****Record loci in each loci group
}

*/

