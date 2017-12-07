#include "Indiv.h"


// This constructor is called only for the individuals in the first generation
Indiv::Indiv( vector < vector<int> >& groupedLoci, double (&deltaS)[3]) {
	
	//cout << "Indiv constructor 1" << "\n";
    
	//Initialize recombination rate
	modifierRecomb[0] = INITIALRECOMBINATION;
	modifierRecomb[1] = INITIALRECOMBINATION;
	
	//Initialize selfing rate 
	modifierSelfing[0] = INITIALSELFING;
	modifierSelfing[1] = INITIALSELFING;
	

	//add initial number of mutations assuming mutation-selection balance
	addInitialMutation(groupedLoci);

	
	//calc indiv fitness, selfing rate and recomb rate, important when calculating population values
	//calcSelfingRate must be calculated BEFORE calcFitness_V3 because indivSelfingRate is used to calc pollen discounted fitness (PDFitness)
	calcRecombRate();
	calcSelfingRate();
	
	//If ignoring a group of loci, then clear the initial mutations from that group
	if(IGNORE_LOCI == 1){
		Cm1[2].clear();
		Cm2[2].clear();
	}
	else if (IGNORE_LOCI == 2){
		Cm1[0].clear();
		Cm1[1].clear();
		
		Cm2[0].clear();
		Cm2[1].clear();
	}

	//calc fitness
	calcFitness_V3(deltaS);
	//cout << "delMut: " << delMut[0] << "," << delMut[1] << "," << delMut[2] << "\n";
	//cout << "Cm1: " << Cm1[0].size() << "," << Cm1[1].size() << "," << Cm1[2].size() << "\n";
	//cout << "Cm2: " << Cm2[0].size() << "," << Cm2[1].size() << "," << Cm2[2].size() << "\n";
}


//Constructor to create individual from info from gametes 1 and 2
Indiv::Indiv( vector < vector<int> > gamete1, vector < vector<int> > gamete2, double recomb1, double recomb2, double selfing1, double selfing2){

	//*****Set mutation list and modifier loci values
	Cm1 = gamete1;
	Cm2 = gamete2;
	
	modifierRecomb[0] = recomb1;
	modifierRecomb[1] = recomb2;
	modifierSelfing[0] = selfing1;
	modifierSelfing[1] = selfing2;
	//*****Set mutation list and modifier loci values

	//calc indiv fitness, selfing rate and recomb rate, important when calculating population values
	//calcSelfingRate must be calculated BEFORE calcFitness_V3 because indivSelfingRate is used to calc pollen discounted fitness (PDFitness)
	calcRecombRate();
	calcSelfingRate();
}


//Add initial number of mutations to indiv based on expectation at mutation-selection balance
void Indiv::addInitialMutation(vector < vector<int> >& groupedLoci){
	
	//initialize Cm1,2 with three empty vectors (each one represent one group of loci)
	for(int i = 0; i < 3; i++){
		Cm1.push_back( vector<int> () );
		Cm2.push_back( vector<int> () );
	}

	//Initialize random number generator
	const gsl_rng_type * T;
	gsl_rng * r;
	T = gsl_rng_default;
	r = gsl_rng_alloc (T);
	gsl_rng_default_seed += time(NULL);
	gsl_rng_set(r, gsl_rng_default_seed);
	
	//These variables record the initial average number of mutations per individual
	int numInitDiploitMut_Del = 0;
	int numInitDiploitMut_Ben = 0;
	if(S > 0){
		numInitDiploitMut_Del = gsl_ran_poisson(r, MUTATIONRATE/(0.5*S) );
	}
	if(S_BEN > 0){
		numInitDiploitMut_Ben = gsl_ran_poisson(r, U_BEN/(0.5*S_BEN) );
	}
			
	//*****Prepare for multinomial sampling
	//Variables to create probability bins
	double	groupProb[3];
	double 	cumulativeProb = 0;
	double	Bin[3];
	
	//record num loci in each group
	double 	numLoci[3];
	
	//generate probability Bins for multinomial sampling
	for(int i = 0; i < 3; i++){
		numLoci[i] = groupedLoci[i].size();
		groupProb[i] = numLoci[i]/L;
		cumulativeProb += groupProb[i];
		Bin[i] = cumulativeProb;
		
		//cout << numLoci[i] << "," << groupProb[i] << "," << cumulativeProb << "," << Bin[i] << "\n";
	}
	//*****Prepare for multinomial sampling
	
	
	
	//countMut counts number of mutations chosen
	int countMut;
	//chosenLoci_Cm record what loci already has mutation in Cm1,2
	vector<int> chosenLoci_Cm1;
	vector<int> chosenLoci_Cm2;
	
	//*****Add new mut to Individual
	//Choose mutations based on multinomial sampling from set of locus X,Y,Z
	//Prevents mutation from being added as duplicates WITHIN a Cm (allow to be same loci between Cm to allow for homozygous mutations)
	
	//set countMut to 0 and clear list of chosenLoci
	//	chosenLoci are always positive integers
	//		never negative even if there are beneficials because it just tracks which loci has gained a mutation (regardless of del or ben) so I won't mutate a loci twice
	//		positive and negative integer is placed into Cm1[] and Cm2[] to indicate type of mutation
	countMut = 0;
	chosenLoci_Cm1.clear();
	chosenLoci_Cm2.clear();

	//Add new mutations until countMut == (numInitDiploitMut_Del+numInitDiploitMut_Ben)
	while( countMut < (numInitDiploitMut_Del+numInitDiploitMut_Ben) ){
		
		//newly chosen mutation
		int tempNewMut;
		//indicate group of newly chosen mut
		int newLocus[3];
		//duplicate loci indicator, initialize as 0
		int duplicateLoci;
		//random chose which Cm to add new mutation
		int whichCm;
		
		
		//do loop to chose new mut AND which Cm to add it to
		//repeat while duplicateLoci == 1 (if current new mut matches already existing new mut)
		do{
			gsl_rng_default_seed += time(NULL)^countMut;
			gsl_rng_set(r, gsl_rng_default_seed);
			//initially assume duplicate is unique
			duplicateLoci = 0;
			//reset new class of new mut
			newLocus[0] = 0; 
			newLocus[1] = 0;
			newLocus[2] = 0;
			
			//Add del mut to any of the three loci groups, but beneficials to only group 2 (constant selection)
			if(countMut < numInitDiploitMut_Del){		
				//uniform random variate
				double UnifRV = gsl_ran_flat(r, 0.0, 1.0);
		
				//choose mut from loci_X if fall between 0 and Bin[0]	
				if( UnifRV <= Bin[0]){
					int tempIndex = floor(gsl_ran_flat(r, -0.5, (numLoci[0]-1)+0.4999) + 0.5);
					tempNewMut = groupedLoci[0][tempIndex];
					newLocus[0] = 1;
				}
				//choose mut from loci_Y if fall between Bin[0] and Bin[1]
				else if ( UnifRV > Bin[0] && UnifRV <= Bin[1]){
					int tempIndex = floor(gsl_ran_flat(r, -0.5, (numLoci[1]-1)+0.4999) + 0.5);
					tempNewMut = groupedLoci[1][tempIndex];
					newLocus[1] = 1;
				}
				//choose mut from loci_Z if fall between Bin[1] and 1
				else{
					int tempIndex = floor(gsl_ran_flat(r, -0.5, (numLoci[2]-1)+0.4999) + 0.5);
					tempNewMut = groupedLoci[2][tempIndex];
					newLocus[2] = 1;
				}
			}
			else{
				int tempIndex = floor(gsl_ran_flat(r, -0.5, (numLoci[2]-1)+0.4999) + 0.5);
				tempNewMut = groupedLoci[2][tempIndex];
				newLocus[2] = 1;
			}
			
			
			//random choose which Cm to check for duplicate
			whichCm = floor(gsl_ran_flat(r, 0.5, 2.4999) + 0.5);
						
			//check through list of chosenLoci_Cm1 or chosenLoci_Cm2 to check if mutation already exists on a specific Cm
			//duplicateLoci == 1 if new chosen mut already exist on Cm1 or Cm2
			if(whichCm == 1){
				vector<int>::iterator itI = chosenLoci_Cm1.begin();
				while( duplicateLoci == 0 && itI < chosenLoci_Cm1.end() ){
					if( tempNewMut == *itI ){
						duplicateLoci = 1;
					}
					itI++;
				}
			}
			else if(whichCm == 2){
				vector<int>::iterator itI = chosenLoci_Cm2.begin();
				while( duplicateLoci == 0 && itI < chosenLoci_Cm2.end() ){
					if( tempNewMut == *itI ){
						duplicateLoci = 1;
					}
					itI++;
				}
			}
		
		}while(duplicateLoci == 1);
		
		//after the above do loop, we have chosen a new mut to be added to a specifc Cm and is NOT a duplicate on that Cm
		
		
		//given new mut is non duplicate within a Cm, add to correct mut class vector (X,Y,Z and del or ben); ben loci multiple by -1
		if(whichCm == 1){
			if(newLocus[0] == 1){
				if(countMut < numInitDiploitMut_Del){
					Cm1[0].push_back(tempNewMut);
				}
				else{
					Cm1[0].push_back(tempNewMut * -1);
				}
			}
			else if(newLocus[1] == 1){
				if(countMut < numInitDiploitMut_Del){
					Cm1[1].push_back(tempNewMut);
				}
				else{
					Cm1[1].push_back(tempNewMut * -1);
				}
			}
			else if(newLocus[2] == 1){
				if(countMut < numInitDiploitMut_Del){
					Cm1[2].push_back(tempNewMut);
				}
				else{
					Cm1[2].push_back(tempNewMut * -1);
				}
			}
			
			//add tempNewMut to vector of chosenLoci_Cm1 (used to check for duplicates)
			chosenLoci_Cm1.push_back(tempNewMut);
		}
		else if (whichCm == 2){
			if(newLocus[0] == 1){
				if(countMut < numInitDiploitMut_Del){
					Cm2[0].push_back(tempNewMut);
				}
				else{
					Cm2[0].push_back(tempNewMut * -1);
				}
			}
			else if(newLocus[1] == 1){
				if(countMut < numInitDiploitMut_Del){
					Cm2[1].push_back(tempNewMut);
				}
				else{
					Cm2[1].push_back(tempNewMut * -1);
				}
			}
			else if(newLocus[2] == 1){
				if(countMut < numInitDiploitMut_Del){
					Cm2[2].push_back(tempNewMut);
				}
				else{
					Cm2[2].push_back(tempNewMut * -1);
				}
			}
			
			//add tempNewMut to vector of chosenLoci_Cm2 (used to check for duplicates)
			chosenLoci_Cm2.push_back(tempNewMut);
		}
		
		//increase count of new mutations
		countMut++;
		
	}
	//*****Add new mut for Individual
	
	//*****Sort mut based on absolute value
	sort(Cm1[0].begin(),Cm1[0].end(),compareAbsolute);
	sort(Cm1[1].begin(),Cm1[1].end(),compareAbsolute);
	sort(Cm1[2].begin(),Cm1[2].end(),compareAbsolute);
	sort(Cm2[0].begin(),Cm2[0].end(),compareAbsolute);
	sort(Cm2[1].begin(),Cm2[1].end(),compareAbsolute);
	sort(Cm2[2].begin(),Cm2[2].end(),compareAbsolute);
		
	//free memory; prevent memory leaks
	gsl_rng_free(r);
}

//Calc indiv selfing rate, ALWAYS called when indiv is created
void Indiv::calcSelfingRate(){
	indivSelfingRate = (modifierSelfing[0] + modifierSelfing[1])/2;
}

//Calc indiv recomb rate, ALWAYS called when indiv is created
void Indiv::calcRecombRate(){
	indivRecombRate = (modifierRecomb[0] + modifierRecomb[1])/2;
}

//calc indivFitness AND indivPDFintess, ALWAYS called when indiv is created
void Indiv::calcFitness_V3( double (&deltaS)[3] ){
	//initialize fitness
	indivFitness = 1;

	//loop through all 3 loci groups
	for(int i = 0; i < 3; i++){
		//size of Cm1 and Cm2
		int sizeCm1 = Cm1[i].size();
		int sizeCm2 = Cm2[i].size();
		//index for Cm1 and Cm2 vector
		int iCm1 = 0;
		int iCm2 = 0;
		//focal mutation in Cm1 and Cm2 being examined
		int mutCm1;
		int mutCm2;
		
		//*****Count number of del and ben loci in Cm1 + Cm2 
		//count number of Del and Ben loci
		delMut[i] = 0;
		benMut[i] = 0;
	
		//Only count ben mut if U_BEN >0
		//	ben mut are negative in value
		if(U_BEN > 0){
			for( vector<int>::iterator itI = Cm1[i].begin(); itI < Cm1[i].end(); itI++ ){
				if( *itI < 0 ){
					benMut[i]++;
				}
			}
			for( vector<int>::iterator itI = Cm2[i].begin(); itI < Cm2[i].end(); itI++ ){
				if( *itI < 0 ){
					benMut[i]++;
				}
			}
		}
		
		//num del loci is total subtract num ben loci
		delMut[i] = sizeCm1 + sizeCm2 - benMut[i];
		//*****Count number of del and ben loci in Cm1 + Cm2 
	
		delHomo[i] = 0;
		benHomo[i] = 0;
		delHetero[i] = 0;
		benHetero[i] = 0;
	
		//*****Find number of homozygous and heterozygous loci (del and ben)
		//loop exit when searched through all of EITHER Cm (where it is no longer possible to have homozygous loci)
		while( iCm1 < sizeCm1 && iCm2 < sizeCm2){
			//mut on Cm1 and Cm2 for comparison
			mutCm1 = Cm1[i][iCm1];
			mutCm2 = Cm2[i][iCm2];
			
			if( abs(mutCm1) < abs(mutCm2) ){
				iCm1++;
			}
			else if (abs(mutCm1) > abs(mutCm2) ){
				iCm2++;
			}
			else if ( abs(mutCm1) == abs(mutCm2) ) {
				//only increase homozygous if sign of mutCm1 and mutCm2 are same, otherwise, move to next mutations
				if(mutCm1 > 0 && mutCm2 > 0){	
					delHomo[i]++;
				}
				else if( mutCm1 < 0 && mutCm2 < 0 ){	
					benHomo[i]++;
				}		
				iCm1++;
				iCm2++;
			}
		}
	
		//calculate num of heterozygous loci
		delHetero[i] = delMut[i] - 2 * delHomo[i];
		benHetero[i] = benMut[i] - 2 * benHomo[i];
		//*****Find number of homozygous and heterozygous loci (del and ben)
		
		//calc cumulative fitness (cumulative product over all i loci groups)
		indivFitness *=	pow( (1 - S * deltaS[i]), delHomo[i] ) *
						pow( (1 - S * H * deltaS[i]), delHetero[i] ) *
						pow( (1 + S_BEN * deltaS[i]), benHomo[i] ) *
						pow( (1 + S_BEN * H_BEN * deltaS[i]), benHetero[i] );
	}
	
	//indivPDFitness depends on selfing rate and intensity of pollen discounting (K)
	indivPDFitness = (1 - K * indivSelfingRate) * indivFitness;
	
	/*
	cout << "Cm:\n";
	PrintVecInt2D(Cm1);
	PrintVecInt2D(Cm2);
	cout << delMut[0] <<","<< delMut[1] <<","<< delMut[2]<<"\n";
	cout << benMut[0] <<","<< benMut[1] <<","<< benMut[2]<<"\n";
	*/
	
	
}




//Takes gamete as reference argument and inserts mutation into it after recombination between Cm1 and Cm2
//First generate list of break points for recombination
//	listBreakPoints is a vector of integers to indicate the break point is directly UPSTREAM of the specified loci[i] and downstream of the loci[i-1]
//	we generate break points in the range [2, L]  since break points are upstream of chosen loci 
//	i.e. break point at locus L POSSIBLE since there is still loci L downstream, break point at locus 1 IMPOSSIBLE since there is no loci upstream of locus 1
//Recombinant genotype is stored in indiv private variable: gamete, gameteModifierRecomb and gameteModifierSelfing
//*****Aug 30 2016, input meanRecombPopul, but actually is not being used; currently just using indivRecombRate to determine numBreakPoints
void Indiv::makeGamete(double& meanRecombPopul, vector< vector<int> >& gamete, double& gameteModifierRecomb, double& gameteModifierSelfing, int (&gameteNeutralAllele)[3]){

	//Initialize RNG	
	const gsl_rng_type * T;
	gsl_rng * r;
	T = gsl_rng_taus;
	r = gsl_rng_alloc (T);
	gsl_rng_default_seed += time(NULL);
	gsl_rng_set (r, gsl_rng_default_seed);
	
	//*****Randomly choose number and list of break points for recombination
	//integer indicating number of break points required
	int numBreakPoints;
	//vector to hold list of randomly chosen break points
	vector<int> listBreakPoints;

	//Poisson distributed number of break points for recombination
	numBreakPoints = gsl_ran_poisson(r, indivRecombRate);
		
	//randomly select numBreakPoints break points and store in listBreakPoints (prevents duplicate break points)
	for(int i = 0; i < numBreakPoints; i++){
		//duplicate indicator
		int duplicateBP;
		//temp break point variable (only store in listBreakPoints if not a duplicate)
		int tempBP;

		//do loop exits only if tempBP is unique within listBreakPoints (loops if duplicateBP gets changed from 0 to 1)
		do{
			//random select break points [2, L]
			tempBP = floor( gsl_ran_flat(r, 1.5, L+0.4999) + 0.5 );
			//assume no duplicate initially
			duplicateBP = 0;
			
			//loop through listBreakPoints to check if tempBP already chosen
			//exit when at end of listBreakPoints OR duplicateLoci == 1 (i.e. tempBP already exist in listBreakPoints)
			vector<int>::iterator itI = listBreakPoints.begin();
			while( duplicateBP == 0 && itI < listBreakPoints.end() ){
				if( tempBP == *itI ){
					duplicateBP = 1;
				}
				itI++;
			}
		}while(duplicateBP == 1);
		
		//add tempBP to listBreakPoints after it had passed duplicate check
		listBreakPoints.push_back(tempBP);

	}
	//sort break points in ascending order
	sort(listBreakPoints.begin(), listBreakPoints.end());
	
	//include 1 and L+1 at ends of listBreakPoints to accomodate algorithm (i.e. algorithm identify segments as listBreakPoints[i] to listBreakPoints[i+1])
	listBreakPoints.insert(listBreakPoints.begin(), 1);
	listBreakPoints.push_back(L+1);
	//*****Randomly choose number and list of break points for recombination
	
	
	//*****Randomly select which Cm will create first segment of gamete
	//update seed
	gsl_rng_default_seed += time(NULL);
	gsl_rng_set (r, gsl_rng_default_seed);
	
	//First segment inherit from Cm1 (2) if focalCm = 0 (1); save first focal Cm in initFocalCm
	int initFocalCm = floor(gsl_ran_flat(r, -0.5, 1.4999) + 0.5);
	int focalCm = initFocalCm;
	//*****Randomly select which Cm will create first segment of gamete
		
		
	//******Create gamete by copy segments of Cm1 and Cm2 to gamete in addition to modifier loci
	//		Gamete info is then stored in private variables: gamete, gameteModifierRecomb, gameteModifierSelfing 
	
	//clear any existing elements in gamete
	gamete.clear();
	
	//initialize gamete with three empty vectors (each one represent one group of loci)
	for(int i = 0; i < 3; i++){
		gamete.push_back( vector<int> () );
	}
	
	vector<int> iCm1(3,0);
	vector<int> iCm2(3,0);
		
	//record number of mutations in each group of Cm1 and Cm2
	vector<int> sizeCm1;
	vector<int> sizeCm2;
	for(int i = 0; i < 3; i++){
		sizeCm1.push_back( Cm1[i].size() );
		sizeCm2.push_back( Cm2[i].size() );
	}
		
	int iBP = 0;
	int segStart;
	int segEnd;
	
	//specify position of selfing and recomb modifier loci
	double modLocusRecomb = floor(L/3.0);
	double modLocusSelfing = floor(2.0*L/3.0);

	//*****Eddie MODDED Mar 29 2016
	//Specify position of neutral loci
	double neutralLocus[3];
	neutralLocus[0] = floor(L/4);
	neutralLocus[1] = floor(L/2);
	neutralLocus[2] = floor(L*3/4);
	//cout << "Loci positions: " << modLocusRecomb << "," << modLocusSelfing << "," << neutralLocus[0] << "," << neutralLocus[1] << "," << neutralLocus[2] << "\n";
	//*****Eddie MODDED Mar 29 2016
	

	
	//copy segments of Cm1 or Cm2 to gamete depending on value at focalCm
	//loop exit when reach last segment of break points ( segEnd >= L+1 )
	do{				
		//update start and end of segment (loop exits before segEnd go out of bounds of listBreakPoints)
		segStart = listBreakPoints[iBP];
		segEnd = listBreakPoints[iBP+1];
		//cout << segStart << "," << segEnd << "\n";
		
		//Copy from Cm1
		if( focalCm == 0){
			//loop through all 3 mut groups
			for(int j = 0; j < 3; j++){
				if( sizeCm1[j] > 0 ){
									
					//go through Cm1[j] and copy mut in Cm1 to gamete belonging to region [segStart, segEnd) (exclusive of segEnd)
					while( abs(Cm1[j][ iCm1[j] ]) < segEnd && iCm1[j] < sizeCm1[j] ){

						if( abs(Cm1[j][ iCm1[j] ]) >= segStart && abs(Cm1[j][ iCm1[j] ]) < segEnd ){
							gamete[j].push_back( Cm1[j][ iCm1[j] ] );
						}
						iCm1[j]++;
					}
				}
			}
			
			//Get modifier alleles if within segment
			if( modLocusRecomb >= segStart && modLocusRecomb < segEnd ){
				gameteModifierRecomb = modifierRecomb[0];
			}
			if( modLocusSelfing >= segStart && modLocusSelfing < segEnd ){
				gameteModifierSelfing = modifierSelfing[0];
			}
			
			//Loop through all three neutral loci
			for(int j = 0; j < 3; j++){
				//Get neutral alleles if within segment
				if( neutralLocus[j] >= segStart && neutralLocus[j] < segEnd ){
					gameteNeutralAllele[j] = 0;
				}
			}
			
			//move to next segment
			focalCm = 1;
			iBP++;
		}
		//Copy from Cm2
		else if( focalCm == 1){	
			//loop through all 3 mut groups
			for(int j = 0; j < 3; j++){
				if( sizeCm2[j] > 0 ){				
					while( abs(Cm2[j][ iCm2[j] ]) < segEnd && iCm2[j] < sizeCm2[j] ){												
						if( abs(Cm2[j][ iCm2[j] ]) >= segStart && abs(Cm2[j][ iCm2[j] ]) < segEnd ){
							gamete[j].push_back( Cm2[j][ iCm2[j] ] );
						}
						iCm2[j]++;						
					}
				}
			}
			//Get modifier alleles if within segment
			if( modLocusRecomb >= segStart && modLocusRecomb < segEnd ){
				gameteModifierRecomb = modifierRecomb[1];
			}
			if( modLocusSelfing >= segStart && modLocusSelfing < segEnd ){
				gameteModifierSelfing = modifierSelfing[1];
			}
			
			//Loop through all three neutral loci
			for(int j = 0; j < 3; j++){
				//Get neutral alleles if within segment
				if( neutralLocus[j] >= segStart && neutralLocus[j] < segEnd ){
					gameteNeutralAllele[j] = 1;
				}
			}
			
			//move to next segment
			focalCm = 0;
			iBP++;
		}
	
	}while( segEnd < (L+1) );
	
	//******Create gamete by copy segments of Cm1 and Cm2 to gamete in addition to modifier loci

	//free memory; prevent memory leaks
	gsl_rng_free(r);
	
}









//*****Eddie MODDED Aug 1 2016
//New constructors and functions for simulations that only want to track neutral loci values
//Default constructor
Indiv::Indiv() {
	
	cout << "Indiv constructor 0" << "\n";
    
	//Initialize recombination rate
	modifierRecomb[0] = INITIALRECOMBINATION;
	modifierRecomb[1] = INITIALRECOMBINATION;
	
	//Initialize selfing rate 
	modifierSelfing[0] = INITIALSELFING;
	modifierSelfing[1] = INITIALSELFING;
	
	//calc selfing rate and recomb rate, important when calculating population values
	calcRecombRate();
	calcSelfingRate();

}



//Constructor to create individual from info from gametes
Indiv::Indiv(double recomb1, double recomb2, double selfing1, double selfing2){


	modifierRecomb[0] = recomb1;
	modifierRecomb[1] = recomb2;
	
	modifierSelfing[0] = selfing1;
	modifierSelfing[1] = selfing2;
	
	//calc selfing rate and recomb rate, important when calculating population values
	calcRecombRate();
	calcSelfingRate();
}


//double (&gameteNeutral)[3] indicates whether Cm1 (0) or Cm2 (1) will be inherited, it does NOT contain the neutral value
void Indiv::makeGamete_New( double& gameteModifierRecomb, double& gameteModifierSelfing, int (&gameteNeutralAllele)[3] ){

	//Initialize RNG	
	const gsl_rng_type * T;
	gsl_rng * r;
	T = gsl_rng_taus;
	r = gsl_rng_alloc (T);
	gsl_rng_default_seed += time(NULL);
	gsl_rng_set (r, gsl_rng_default_seed);
	
	
	//*****Randomly choose number and list of break points for recombination
	//integer indicating number of break points required
	int numBreakPoints;
	//vector to hold list of randomly chosen break points
	vector<int> listBreakPoints;

	//Poisson distributed number of break points for recombination
	numBreakPoints = gsl_ran_poisson(r, indivRecombRate);
		
	//randomly select numBreakPoints break points and store in listBreakPoints (prevents duplicate break points)
	for(int i = 0; i < numBreakPoints; i++){
		//duplicate indicator
		int duplicateBP;
		//temp break point variable (only store in listBreakPoints if not a duplicate)
		int tempBP;

		//do loop exits only if tempBP is unique within listBreakPoints (loops if duplicateBP gets changed from 0 to 1)
		do{
			//random select break points [2, L]
			tempBP = floor( gsl_ran_flat(r, 1.5, L+0.4999) + 0.5 );
			//assume no duplicate initially
			duplicateBP = 0;
			
			//loop through listBreakPoints to check if tempBP already chosen
			//exit when at end of listBreakPoints OR duplicateLoci == 1 (i.e. tempBP already exist in listBreakPoints)
			vector<int>::iterator itI = listBreakPoints.begin();
			while( duplicateBP == 0 && itI < listBreakPoints.end() ){
				if( tempBP == *itI ){
					duplicateBP = 1;
				}
				itI++;
			}
		}while(duplicateBP == 1);
		
		//add tempBP to listBreakPoints after it had passed duplicate check
		listBreakPoints.push_back(tempBP);

	}
	//sort break points in ascending order
	sort(listBreakPoints.begin(), listBreakPoints.end());
	
	
	//include 1 and L+1 at ends of listBreakPoints to accomodate algorithm (i.e. algorithm identify segments as listBreakPoints[i] to listBreakPoints[i+1])
	listBreakPoints.insert(listBreakPoints.begin(), 1);
	listBreakPoints.push_back(L+1);
	//*****Randomly choose number and list of break points for recombination
	
	
	//*****Randomly select which Cm will create first segment of gamete
	//update seed
	gsl_rng_default_seed += time(NULL);
	gsl_rng_set (r, gsl_rng_default_seed);
	
	//First segment inherit from Cm1 (2) if focalCm = 0 (1); save first focal Cm in initFocalCm
	int initFocalCm = floor(gsl_ran_flat(r, -0.5, 1.4999) + 0.5);
	int focalCm = initFocalCm;
	//*****Randomly select which Cm will create first segment of gamete
		
		
	//******Create gamete by copy segments of Cm1 and Cm2 to gamete in addition to modifier loci
	//		Gamete info is then stored in variables: gameteModifierRecomb, gameteModifierSelfing 
	int iBP = 0;
	int segStart;
	int segEnd;
	
	//specify position of selfing and recomb modifier loci
	double modLocusRecomb = floor(L/3.0);
	double modLocusSelfing = floor(2.0*L/3.0);

	//Specify position of neutral loci
	double neutralLocus[3];
	neutralLocus[0] = floor(L/4);
	neutralLocus[1] = floor(L/2);
	neutralLocus[2] = floor(L*3/4);
	//cout << "Loci positions: " << modLocusRecomb << "," << modLocusSelfing << "," << neutralLocus[0] << "," << neutralLocus[1] << "," << neutralLocus[2] << "\n";
	

	
	//copy segments of Cm1 or Cm2 to gamete depending on value at focalCm
	//loop exit when reach last segment of break points ( segEnd >= L+1 )
	do{				
		//update start and end of segment (loop exits before segEnd go out of bounds of listBreakPoints)
		segStart = listBreakPoints[iBP];
		segEnd = listBreakPoints[iBP+1];
		
		//Copy from Cm1
		if( focalCm == 0){
			
			//Get modifier alleles if within segment
			if( modLocusRecomb >= segStart && modLocusRecomb < segEnd ){
				gameteModifierRecomb = modifierRecomb[0];
			}
			if( modLocusSelfing >= segStart && modLocusSelfing < segEnd ){
				gameteModifierSelfing = modifierSelfing[0];
			}
			
			//Loop through all three neutral loci
			for(int j = 0; j < 3; j++){
				//Get neutral alleles if within segment
				if( neutralLocus[j] >= segStart && neutralLocus[j] < segEnd ){
					gameteNeutralAllele[j] = 0;
				}
			}
			
			//move to next segment
			focalCm = 1;
			iBP++;
		}
		//Copy from Cm2
		else if( focalCm == 1){	
			
			//Get modifier alleles if within segment
			if( modLocusRecomb >= segStart && modLocusRecomb < segEnd ){
				gameteModifierRecomb = modifierRecomb[1];
			}
			if( modLocusSelfing >= segStart && modLocusSelfing < segEnd ){
				gameteModifierSelfing = modifierSelfing[1];
			}
			
			//Loop through all three neutral loci
			for(int j = 0; j < 3; j++){
				//Get neutral alleles if within segment
				if( neutralLocus[j] >= segStart && neutralLocus[j] < segEnd ){
					gameteNeutralAllele[j] = 1;
				}
			}
			
			//move to next segment
			focalCm = 0;
			iBP++;
		}
	
	}while( segEnd < (L+1) );
	
	//******Create gamete by copy segments of Cm1 and Cm2 to gamete in addition to modifier loci


	//free memory; prevent memory leaks
	gsl_rng_free(r);
	
}

//*****Eddie MODDED Aug 1 2016




//calc how many homozygoes and heterozygous sites are in the individual
void Indiv::calcGenotype_V3(){

	//loop through all 3 loci groups
	for(int i = 0; i < 3; i++){
		//size of Cm1 and Cm2
		int sizeCm1 = Cm1[i].size();
		int sizeCm2 = Cm2[i].size();
		//index for Cm1 and Cm2 vector
		int iCm1 = 0;
		int iCm2 = 0;
		//focal mutation in Cm1 and Cm2 being examined
		int mutCm1;
		int mutCm2;
		
		//*****Count number of del and ben loci in Cm1 + Cm2 
		//count number of Del and Ben loci
		delMut[i] = 0;
		benMut[i] = 0;
	
		//Only count ben mut if U_BEN >0
		if(U_BEN > 0){
			for( vector<int>::iterator itI = Cm1[i].begin(); itI < Cm1[i].end(); itI++ ){
				if( *itI < 0 ){
					benMut[i]++;
				}
			}
			for( vector<int>::iterator itI = Cm2[i].begin(); itI < Cm2[i].end(); itI++ ){
				if( *itI < 0 ){
					benMut[i]++;
				}
			}
		}
		
		//num del loci is total subtract num ben loci
		delMut[i] = sizeCm1 + sizeCm2 - benMut[i];
		//*****Count number of del and ben loci in Cm1 + Cm2 
	
		delHomo[i] = 0;
		benHomo[i] = 0;
		delHetero[i] = 0;
		benHetero[i] = 0;
	
		//*****Find number of homozygous and heterozygous loci (del and ben)
		//loop exit when searched through all of EITHER Cm (where it is no longer possible to have homozygous loci)
		while( iCm1 < sizeCm1 && iCm2 < sizeCm2){
			//mut on Cm1 and Cm2 for comparison
			mutCm1 = Cm1[i][iCm1];
			mutCm2 = Cm2[i][iCm2];
			
			if( abs(mutCm1) < abs(mutCm2) ){
				iCm1++;
			}
			else if (abs(mutCm1) > abs(mutCm2) ){
				iCm2++;
			}
			else if ( abs(mutCm1) == abs(mutCm2) ) {
				//only increase homozygous if sign of mutCm1 and mutCm2 are same, otherwise, move to next mutations
				if(mutCm1 > 0 && mutCm2 > 0){	
					delHomo[i]++;
				}
				else if( mutCm1 < 0 && mutCm2 < 0 ){	
					benHomo[i]++;
				}		
				iCm1++;
				iCm2++;
			}
		}
	
		//calculate num of heterozygous loci
		delHetero[i] = delMut[i] - 2 * delHomo[i];
		benHetero[i] = benMut[i] - 2 * benHomo[i];
		//*****Find number of homozygous and heterozygous loci (del and ben)
		
	}
}


