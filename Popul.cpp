#include "Popul.h"



//Constructor: When an object of type "Popul" is created, it automatically creates its first generation
Popul::Popul( vector < vector<int> >& groupedLoci, double (&deltaS)[3] ) {

	//Construct individuals and save in thePop
	for(int i=0; i<POPULATION; i++) {
		gsl_rng_default_seed += time(NULL)^i;
		
		//create indiv newIndiv (automaticcally calc its fitness given deltaS array
		Indiv newIndiv(groupedLoci, deltaS);
		
		//save newIndiv in thePop (vector of indiv)
		thePop.push_back(newIndiv);
	}
	
	//Update population measurements
	findHighestFitness();	
	calcRecombRate();
	
	//Initialize fixations as 0
	for(int g = 0; g < 3; g++){	
		fixBen[g] = 0;
		fixDel[g] = 0;
	}
	
	//*****Eddie MODDED Aug, 30, 2016
	//initialize neutralGameteIndex with POPULATION,3,2,2 elements in dimension 1,2,3,4, respectively
	neutralGameteIndex.resize(POPULATION);
	for(int n = 0; n < POPULATION; n++){
		neutralGameteIndex[n].resize(3);
		
		for(int l = 0; l < 3; l++){
			neutralGameteIndex[n][l].resize(2);
			
			for(int a = 0; a < 2; a++){
				neutralGameteIndex[n][l][a].resize(2);
				
				//initialize the parental indices as -1; will give error unless the value is changed to 0 or 1
				neutralGameteIndex[n][l][a][0]=-1;
				neutralGameteIndex[n][l][a][1]=-1;
			}
		}
	}
	//*****Eddie MODDED Aug, 30, 2016
}

//Create new generation of individuals by calling makeOffspring
//thePop replaced by individuals in thePopOffspring and population measurements are updated
//Requires groupedLoci since involves adding new mutations (requires groupedLoci)
void Popul::makeNewGeneration( vector < vector<int> >& groupedLoci ){
	
	//Initialize summary statistic of parent ID to 0 each generation	
	ID_P1 = 0;
	ID_P2 = 0;
	
	//create POPULATION number of offspring stored into thePopOffspring
	for(int i = 0; i < POPULATION; i++) {
		makeOffspring( groupedLoci, i );
	}
	
	//Calc average parental ID
	ID_P1/=POPULATION;
	ID_P2/=POPULATION;
	
	//Update thePop as thePopOffspring which creates the next generation
	thePop = thePopOffspring;
		
	//Clear thePopOffspring, extremely important when making offspring for the generation after this
	thePopOffspring.clear();
		
	//Update population measurements
	calcRecombRate();
}


//Assumes indivFitness and indivPDFitness already calculated
void Popul::findHighestFitness(){
	
	//Initialize highestFitness and highestPDFitness with indiv 0
	populHighestFitness = thePop[0].indivFitness;
	populHighestPDFitness = thePop[0].indivPDFitness;
	
	//loop through i = 1 to POPULATION-1
	for(int i = 1; i < POPULATION; i++){
		
		//if current highestFitness < indivFitness of indiv i, update highestFitness to that of indiv i
		if( populHighestFitness < thePop[i].indivFitness){
			populHighestFitness = thePop[i].indivFitness;
		}
		if( populHighestPDFitness < thePop[i].indivPDFitness){
			populHighestPDFitness = thePop[i].indivPDFitness;
		}	
	}
}

//Calc mean and var individual fitness
void Popul::calcPopulFitness(){
	
	//Initialize populFitness = 0
	populFitness[0] = 0;
	populFitness[1] = 0;
	
	long double sumSquares = 0;
	
	//Get values from all ind
	for(int i = 0; i < POPULATION; i++){
		populFitness[0] += thePop[i].indivFitness;
		sumSquares += powl( thePop[i].indivFitness, 2.0 );
	}
	
	//calc mean and var
	populFitness[0] /= POPULATION;
	populFitness[1] = sumSquares/POPULATION - powl( populFitness[0], 2.0 );
}

//Assumes indivRecombRate already calculated
void Popul::calcRecombRate(){

	//Initialize populRecombRate = 0
	populRecombRate[0] = 0;
	populRecombRate[1] = 0;
	
	double sumSquares = 0;
	
	//Get values from all ind
	for(int i = 0; i < POPULATION; i++){
		populRecombRate[0] += thePop[i].indivRecombRate;
		sumSquares += pow( thePop[i].indivRecombRate, 2.0 );
	}
	
	//calc mean and var
	populRecombRate[0] /= POPULATION;
	populRecombRate[1] = sumSquares/POPULATION - pow( populRecombRate[0], 2.0 );
}

//Assumes indivRecombRate already calculated
void Popul::calcSelfingRate(){

	//Initialize populSelfingRate = 0
	populSelfingRate[0] = 0;
	populSelfingRate[1] = 0;
	
	double sumSquares = 0;
	
	//Get values from all ind
	for(int i = 0; i < POPULATION; i++){
		populSelfingRate[0] += thePop[i].indivSelfingRate;
		sumSquares += pow( thePop[i].indivSelfingRate, 2.0 );
	}
	
	//calc mean and var
	populSelfingRate[0] /= POPULATION;
	populSelfingRate[1] = sumSquares/POPULATION - pow( populSelfingRate[0], 2.0 );
}


//select index of first parent for producing offspring
//return index of parent that pass selection
int Popul::pickFirstParent(){
	
	//index of indiv in thePop that passes selection
	int indexParent;
	//random uniform variate
	long double requiredFitness;
	
	//Initialize RNG	
	const gsl_rng_type * T;
	gsl_rng * r;
	T = gsl_rng_taus;
	r = gsl_rng_alloc (T);
	gsl_rng_default_seed += time(NULL);
	gsl_rng_set (r, gsl_rng_default_seed);
		
	//continue to pick parent until parent passes fitness check
	do{
		//random pick index from 0 to POPULATION-1 (for population with size N, it is indexed from 0 to N-1)
		indexParent = floor(gsl_ran_flat(r, -0.5, (POPULATION-1)+0.4999) + 0.5);
				
		//randomization (not really needed)
		gsl_rng_default_seed += time(NULL)^indexParent;
		gsl_rng_set (r, gsl_rng_default_seed);
		
		//random uniform variate
		requiredFitness = gsl_ran_flat(r, 0, 1);
				
	}while( requiredFitness > (thePop[indexParent].indivFitness/populHighestFitness) );

	//free memory; prevent memory leaks
	gsl_rng_free(r);

	//return index of indiv that passed fitness check
	return indexParent;
}


//select index of parent for outcrossing, utilizes pollen discounted fitness (indivPDFitness)
//return index of parent that pass selection based on pollen discounted fitness
int Popul::pickOutcrossParent( int indexParent1 ){
	
	//index of indiv in thePop that passes selection
	int indexParent2;
	//random uniform variate
	long double requiredFitness;
	
	//Initialize RNG	
	const gsl_rng_type * T;
	gsl_rng * r;
	T = gsl_rng_taus;
	r = gsl_rng_alloc (T);
	gsl_rng_default_seed += time(NULL);
	gsl_rng_set (r, gsl_rng_default_seed);

	//continue to pick parent until parent passes pollen discounted fitness check
	do{
		
		//random pick index from 0 to POPULATION-1 that is != indexParent1
		do{
			indexParent2 = floor(gsl_ran_flat(r, -0.5, (POPULATION-1)+0.4999) + 0.5);
		}while(indexParent2 == indexParent1);
		
		//randomization (not really needed)
		gsl_rng_default_seed += time(NULL)^indexParent2;
		gsl_rng_set (r, gsl_rng_default_seed);
		
		//random uniform variate
		requiredFitness = gsl_ran_flat(r, 0, 1);
		
	}while( requiredFitness > (thePop[indexParent2].indivPDFitness/populHighestPDFitness) );

	//free memory; prevent memory leaks
	gsl_rng_free(r);

	//return index of indiv that passed pollen discounted fitness check
	return indexParent2;
}


//Create a SINGLE offspring individual
//1. Select parents based on fitness and polen discounted fitness if outcrossing
//2. Parents create gametes with recombination rate of populRecombRate[0]
//3. Add new mutations to gametes
//4. Mutate modifier loci if required
//5. Create Indiv Offspring given gamete1, gamete2 and modifier loci values
//6. Store Offspring in thePopOffspring vector
void Popul::makeOffspring( vector < vector<int> >& groupedLoci, int n ){
	
	//Initialize RNG	
	const gsl_rng_type * T;
	gsl_rng * r;
	T = gsl_rng_taus;
	r = gsl_rng_alloc (T);
	gsl_rng_default_seed += time(NULL);
	gsl_rng_set (r, gsl_rng_default_seed);
	
	//index of parent 1 and parent 2 within thePop
	int indexP1, indexP2;
	
	//*****BEGIN Selection: select parents based on their individual fitness relative to the population highest fitness
	//P1 outcrosses if selfing rate < Uniform random variable [0,1], P2 chosen based on pollen discount fitness
	//P1 selfs if selfing rate >= Uniform random variable [0,1], P2 = P1; (i.e. always selfes if selfing rate == 1)
	
	//select parent 1 based on indivFitness
	indexP1 = pickFirstParent();

	//decide if P1 will self or outcross
	double UnifRV = gsl_ran_flat(r, 0, 1);
	
	//Initialize OffspringID
	int OffspringID = 0;
	
	//P1 outcross with P2 != P1 (P2 selected using pollen discounted fitness)
	if( UnifRV > thePop[indexP1].indivSelfingRate ) {
		indexP2 = pickOutcrossParent(indexP1);
		OffspringID = 1;
	}
	//P1 self, so P2 = P1
	else {
		indexP2 = indexP1;
		OffspringID = 0;
	}
	
	//Add up parent1 and parent2 IDs
	ID_P1 += thePop[indexP1].ID;
	ID_P2 += thePop[indexP2].ID;
	//*****END Selection: select parents based on their individual fitness relative to population highest fitness
	
	//*****BEGIN Reproduction: Create gamete for parent 1,2; individuals are thePop[indexP1] and thePop[indexP2]
	//vector to store gamete of parent 1,2 (NO NEED to initialize with empty vectors, indiv makeGamete function will perform that action)
	vector< vector<int> > gamete1;
	vector< vector<int> > gamete2;
	
	double modifierRecomb_gamete1;
	double modifierSelfing_gamete1;
	double modifierRecomb_gamete2;
	double modifierSelfing_gamete2;

	//Initialize neutral loci values with 0
	int neutral_P1Allele[3]={};
	int neutral_P2Allele[3]={};
	
	//makeGamete for parent 1,2 and store info into appropriate variables
	thePop[indexP1].makeGamete( populRecombRate[0], gamete1, modifierRecomb_gamete1, modifierSelfing_gamete1, neutral_P1Allele );
	thePop[indexP2].makeGamete( populRecombRate[0], gamete2, modifierRecomb_gamete2, modifierSelfing_gamete2, neutral_P2Allele );
	//*****Create gamete for parent 1,2
	
	//*****Add new mutations to gamete1,2
	//IGNORE mutations to fluctuating or constant loci depending on value of IGNORE_LOCI
	
	if(IGNORE_LOCI == 0){
		addMutation_V4( gamete1, gamete2, groupedLoci, 0);
	}
	else if (IGNORE_LOCI == 1){
		addMutation_V4_IGNORE1( gamete1, gamete2, groupedLoci, 0);
	}
	else if (IGNORE_LOCI == 2){
		addMutation_V4_IGNORE2( gamete1, gamete2, groupedLoci, 0);
	}

	//*****Add new mutations to gamete1,2	
	
	
	//*****Mutate recombination modifier loci if required
	if( RECOMBMUTATIONRATE > 0){
	
		gsl_rng_default_seed += time(NULL);
		gsl_rng_set (r, gsl_rng_default_seed);
		
		//Uniform random variable to determine if modifier locus will mutate
		double UnifRV;
		
		UnifRV = gsl_ran_flat(r, 0, 1);
		//mutate modifierRecomb_gamete1 if UnifRV < RECOMBMUTATIONRATE (i.e. always mutate when RECOMBMUTATIONRATE = 1)
		if( UnifRV < RECOMBMUTATIONRATE ){
	
			//recomb modifier change by a Gaussian distributed random variable with mean 0 and variance SIGMAR
			double delta = gsl_ran_gaussian(r, SIGMAR);
		
			//mutate modifierRecomb_gamete1
			modifierRecomb_gamete1 += delta;
			
			//modifierRecomb_gamete1 must be bounded by 0 and 1
			if( modifierRecomb_gamete1 < 0 ){
				modifierRecomb_gamete1 = 0;
			}
			else if( modifierRecomb_gamete1 > 1 ){
				modifierRecomb_gamete1 = 1;
			}
		}
		
		UnifRV = gsl_ran_flat(r, 0, 1);
		//mutate modifierRecomb_gamete2 if UnifRV < RECOMBMUTATIONRATE (i.e. always mutate when RECOMBMUTATIONRATE = 1)
		if( UnifRV < RECOMBMUTATIONRATE ){
	
			//recomb modifier change by a Gaussian distributed random variable with mean 0 and variance SIGMAR
			double delta = gsl_ran_gaussian(r, SIGMAR);
		
			//mutate modifierRecomb_gamete1
			modifierRecomb_gamete2 += delta;
			
			//modifierRecomb_gamete1 must be bounded by 0 and 1
			if( modifierRecomb_gamete2 < 0 ){
				modifierRecomb_gamete2 = 0;
			}
			else if( modifierRecomb_gamete2 > 1 ){
				modifierRecomb_gamete2 = 1;
			}
		}
	
	}
	//*****Mutate recombination modifier loci if required
	
	//*****Mutate selfing modifier loci if required
	if( SELFINGMUTATIONRATE > 0){
	
		gsl_rng_default_seed += time(NULL);
		gsl_rng_set (r, gsl_rng_default_seed);
		
		//Uniform random variable to determine if modifier locus will mutate
		double UnifRV;
		
		UnifRV = gsl_ran_flat(r, 0, 1);
		//mutate modifierRecomb_gamete1 if UnifRV < SELFINGMUTATIONRATE (i.e. always mutate when SELFINGMUTATIONRATE = 1)
		if( UnifRV < SELFINGMUTATIONRATE ){
	
			//recomb modifier change by a Gaussian distributed random variable with mean 0 and variance SIGMAR
			double delta = gsl_ran_gaussian(r, SIGMAS);
		
			//mutate modifierRecomb_gamete1
			modifierSelfing_gamete1 += delta;
			
			//modifierRecomb_gamete1 must be bounded by 0 and 1
			if( modifierSelfing_gamete1 < 0 ){
				modifierSelfing_gamete1 = 0;
			}
			else if( modifierSelfing_gamete1 > 1 ){
				modifierSelfing_gamete1 = 1;
			}
		}
		
		UnifRV = gsl_ran_flat(r, 0, 1);
		//mutate modifierRecomb_gamete1 if UnifRV < SELFINGMUTATIONRATE (i.e. always mutate when SELFINGMUTATIONRATE = 1)
		if( UnifRV < SELFINGMUTATIONRATE ){
	
			//recomb modifier change by a Gaussian distributed random variable with mean 0 and variance SIGMAR
			double delta = gsl_ran_gaussian(r, SIGMAS);
		
			//mutate modifierRecomb_gamete1
			modifierSelfing_gamete2 += delta;
			
			//modifierRecomb_gamete1 must be bounded by 0 and 1
			if( modifierSelfing_gamete2 < 0 ){
				modifierSelfing_gamete2 = 0;
			}
			else if( modifierSelfing_gamete2 > 1 ){
				modifierSelfing_gamete2 = 1;
			}
		}
	}
	//*****Mutate selfing modifier loci if required
	
	
	//*****Create offspring using constructor of Indiv
	Indiv Offspring( gamete1, gamete2, modifierRecomb_gamete1, modifierRecomb_gamete2, modifierSelfing_gamete1, modifierSelfing_gamete2);
	//*****Create offspring using constructor of Indiv	
	
	//*****Nov18
	//Manually assign ID 
	Offspring.ID = OffspringID;
	//*****Nov18	
	
	//*****Store Offspring in thePopOffspring
	//Assumes thePopOffspring has been cleared in the previous generation
	thePopOffspring.push_back( Offspring );
	//*****Store Offspring in thePopOffspring
	
	if(Q_NEUTRAL_LOCI == 1){	
		//*****Record neutral indices
		//For offspring n's allele 0 for all three neutral loci
		neutralGameteIndex[n][0][0][0]=indexP1;
		neutralGameteIndex[n][0][0][1]=neutral_P1Allele[0];
		neutralGameteIndex[n][1][0][0]=indexP1;
		neutralGameteIndex[n][1][0][1]=neutral_P1Allele[1];
		neutralGameteIndex[n][2][0][0]=indexP1;
		neutralGameteIndex[n][2][0][1]=neutral_P1Allele[2];
	
		//For offspring n's allele 1 for all three neutral loci
		neutralGameteIndex[n][0][1][0]=indexP2;
		neutralGameteIndex[n][0][1][1]=neutral_P2Allele[0];
		neutralGameteIndex[n][1][1][0]=indexP2;
		neutralGameteIndex[n][1][1][1]=neutral_P2Allele[1];
		neutralGameteIndex[n][2][1][0]=indexP2;
		neutralGameteIndex[n][2][1][1]=neutral_P2Allele[2];

		/*Display neutralGameteIndex info
		cout << "Print neutralGameteIndex:\n";
		cout << neutralGameteIndex[n][0][0][0] << "," << neutralGameteIndex[n][0][0][1] << "\n";
		cout << neutralGameteIndex[n][1][0][0] << "," << neutralGameteIndex[n][1][0][1] << "\n";
		cout << neutralGameteIndex[n][2][0][0] << "," << neutralGameteIndex[n][2][0][1] << "\n";
		cout << neutralGameteIndex[n][0][1][0] << "," << neutralGameteIndex[n][0][1][1] << "\n";
		cout << neutralGameteIndex[n][1][1][0] << "," << neutralGameteIndex[n][1][1][1] << "\n";
		cout << neutralGameteIndex[n][2][1][0] << "," << neutralGameteIndex[n][2][1][1] << "\n";
		*/
		//*****Record neutral indices	
	}
	//*****END Reproduction
	
	//free memory; prevent memory leaks
	gsl_rng_free(r);
}



//add vMut into vGamete in ascending order
//	assumes vMut and vGamete consists of three vectors representing mut in group 0,1,2
//	assumes no duplicate values within vMut and vGamete (but allow duplicate among vMut and vGamete)
//	assumes vMut[i] and vGamete[i] already sorted in ascending order
void Popul::addMutToGamete(vector< vector<int> >& vMut, vector< vector<int> >& vGamete){
	
	//loop through the 3 vectors within vMut and vGanete
	for(int i = 0; i < 3; i++){
			
		//if vGamete[i] empty (no mutations), insert vMut[i] to vGamete[i]
		if (vGamete[i].size() == 0){
			vGamete[i] = vMut[i];
		}
		//else vGamete[i] non-empty, need to add new mut and keep ascending order
		else{
		
			//index for gamete vector and mutation vector set to 0, will keep track of focal element of vectors 
			int iG = 0;
			int iM = 0;
			//continue adding mut until iM reaches end of vMut (i.e. all mut added)
			while( iM < vMut[i].size() ){
				
				if( abs(vMut[i][iM]) < abs(vGamete[i][iG]) ){
					//add new mut before vGamete[iG], index = vGamete.begin() + iG; .insert() adds to position BEFORE the given element
					vGamete[i].insert( vGamete[i].begin() + iG, vMut[i][iM] );			
					iM++;		
				}
				//if NewMut locus larger than focal locus at Gamete1, continue forward on vGamete
				//	if at last element of vGamete, add rest of vMut to end of vGamete (since rest of vMut larger than last locus of vGamete, given that both were initially sorted ascending)
				else if( abs(vMut[i][iM]) > abs(vGamete[i][iG]) ){
				
					//iG == (vGamete[i].size() - 1), we have reached last element in vGamete
					//	add rest of vMut after the last element of vGamete[i] == vGamete[i].begin() + iG +1 == vGamete[i].begin() + (vGamete[i].size() - 1) +1 (this was Ahmad's code)
					//	more briefly add before one-after-the-last-element == vGamete[i].end() 
					if( iG == (vGamete[i].size() - 1) ){
						vGamete[i].insert( vGamete[i].end(), vMut[i].begin() + iM, vMut[i].end() );
					
						//set iM = vMut[i].size() since all mut added and can exit loop
						iM = vMut[i].size();
					}			
				}
				//vMut[i] locus == locus at vGamete[i], do nothing if same sign, change sign if different sign
				else{
					if ( (vMut[i][iM] * vGamete[i][iG]) < 0 ){
						vGamete[i][iG] *= -1;
					}
					iM++;			
				}
			
				//each of the statements above need iG++ but only until it indexes the last element (vGamete[i].size() - 1)
				if(iG < (vGamete[i].size() - 1)){
					iG++;
				}
			}
		}
	}
}


//Assume all Cm are in ascending order
//Remove all fixed mutations from all 3 mut groups
void Popul::removeFixations(){	
	//*****Find list of fixed mutations, record in listFixations
	//Initialize list of fixed mutations as Cm1 of ind 0
	vector< vector<int> > listFixations;
	listFixations = thePop[0].Cm1;
	
	//total number of fixed mutations	
	int listFixationSize = listFixations[0].size() + listFixations[1].size() + listFixations[2].size();
		
	//compare all Cm1 to listFixations, remove mut in listFixation that do not occur in Cm1 of any individual
	for (int n = 0; n < POPULATION; n++){
		
		//update listFixationSize at each iteration
		listFixationSize = listFixations[0].size() + listFixations[1].size() + listFixations[2].size();
	
		//break out of for loop if listFixations empty
		if( listFixationSize == 0 ){
			break;
		}
	
		//listTemp REFERENCES to individual n's Cm1 for comparing to listFixation
		vector< vector<int> >& listTemp = thePop[n].Cm1;
				
		//loop through all three mutation groups
		for(int g = 0; g < 3; g++){
			
			int iFix = 0;
			int iTemp = 0;
				
			//while loop will compare listFixation and listTemp and remove all mut in listFixation not found in listTemp
			//Exits when reach end of listFixation or end of listTemp
			while( iFix < listFixations[g].size() && iTemp < listTemp[g].size() ){
				
				//element of listFixation or listTemp currently being compared				
				int eFix = listFixations[g][iFix];
				int eTemp = listTemp[g][iTemp];
			
				//If perfect match (sign and value), move to next element in listFixation AND listTemp; retain this element in listFixations
				//	otherwise compare absolute values to see how to proceed with comparisons
				if ( eFix == eTemp ){
					iFix++;
					iTemp++;
				}
				//else if only match in abs value and not sign, remve eFix from listFixation (there will be no other match)
				//NO need to increase iFix but need to move to next element of listTemp
				else if ( abs(eFix) == abs(eTemp) ){
					listFixations[g].erase( listFixations[g].begin() + iFix );
					iTemp++;
				}
				//else if eFix < eTemp, it is impossible for eFix to match listTemp since all elements larger than eFix (assuming all Cm sorted ascending), remove eFix from listFixation; 
				//NO need to increase iFix, since listFixation shortens by one element
				else if( abs(eFix) < abs(eTemp) ){
					listFixations[g].erase( listFixations[g].begin() + iFix );
				}
				//else if eFix > eTemp, still possible to find match for eFix in listTemp, move to next element in listTemp
				else if ( abs(eFix) > abs(eTemp) ){
					iTemp++;
				}
				
			}
			
			//all remaining mut in listFixation will NOT find match if gone through all of listTemp
			//remove all elements in listFixation after while loop ends under violating iTemp < listTemp[g].size()
			//if actually ended by violating iFix < listFixations[g].size(), no mutations are removed since listFixations[g].begin() + iFix already at end of vector
			listFixations[g].erase( listFixations[g].begin() + iFix, listFixations[g].end() );
		}
	}
	
	//compare all Cm2 to listFixations, remove mut in listFixation that do not occur in Cm2 of any individual
	for (int n = 0; n < POPULATION; n++){
	
		//update listFixationSize at each iteration
		listFixationSize = listFixations[0].size() + listFixations[1].size() + listFixations[2].size();
	
		//break out of for loop if listFixations completely empty
		if( listFixationSize == 0 ){
			break;
		}
	
		//listTemp REFERENCES to a focal Cm1 for comparing to listFixation
		vector< vector<int> >& listTemp = thePop[n].Cm2;
	
		//loop through all three mutation groups
		for(int g = 0; g < 3; g++){
		
			int iFix = 0;
			int iTemp = 0;
	
			//while loop will compare listFixation and listTemp and remove all mut in listFixation not found in listTemp
			while( iFix < listFixations[g].size() && iTemp < listTemp[g].size() ){
		
				int eFix = listFixations[g][iFix];
				int eTemp = listTemp[g][iTemp];
				
				//Explanation above when examining Cm1
				if ( eFix == eTemp ){
					iFix++;
					iTemp++;
				}
				else if ( abs(eFix) == abs(eTemp) ){
					listFixations[g].erase( listFixations[g].begin() + iFix );
					iTemp++;
				}
				else if( abs(eFix) < abs(eTemp) ){
					listFixations[g].erase( listFixations[g].begin() + iFix );
				}
				else if ( abs(eFix) > abs(eTemp) ){
					iTemp++;
				}
			}

			listFixations[g].erase( listFixations[g].begin() + iFix, listFixations[g].end() );
		}
	}
	//*****Find list of fixed mutations, record in listFixations

	//Final update for listFixationSize
	listFixationSize = listFixations[0].size() + listFixations[1].size() + listFixations[2].size();
		
	//*****Eddie MODDED Mar 7 2016
	//Checking if fixed mutations indicated in listFixations actually exist in ALL Cm of ALL ind
	if(CHECK_FIXATION == 1 & listFixationSize>0){	
		checkFixations(listFixations);
	}
	//*****Eddie MODDED Mar 7 2016
	
	//*****Record and remove fixations if listFixationSize > 0
	vector<int> newFixBen(3,0);
	if( listFixationSize > 0){
	
		//*****Record number of fixations in each mutation group
	
		//count beneficial mutation if U_BEN > 0
		if(U_BEN > 0){
			for(int g = 0; g < 3; g++){
				for( vector<int>::iterator itI = listFixations[g].begin(); itI < listFixations[g].end(); itI++ ){
					if( *itI < 0 ){
						fixBen[g]++;
						newFixBen[g]++;
					}
				}
			}
		}
	
		//record number of fixed deleterious mut
		fixDel[0] += listFixations[0].size() - newFixBen[0];
		fixDel[1] += listFixations[1].size() - newFixBen[1];
		fixDel[2] += listFixations[2].size() - newFixBen[2];
	
		//*****Record number of fixations in each mutation group
		
		//*****Remove fixed mutation from all individuals and Cm
		//loop through all ind in population
		for (int n = 0; n < POPULATION; n++){
			
			//reference tempCm1 and tempCm2 as Cm of ind n
			vector< vector<int> >& tempCm1 = thePop[n].Cm1;
			vector< vector<int> >& tempCm2 = thePop[n].Cm2;
		
			//loop through all 3 mutation groups
			for(int g = 0; g < 3; g++){
			
				//remove fixations from Cm1
				int iFix = 0;
				int iTemp = 0;
				while ( iFix < listFixations[g].size() ){
					//remove tempCm1[g][iTemp] if it is equal listFixation[g][iFix]
					//then move to next element in listFixation (iFix++) and tempCm1 (no need to increase iTemp after erasing an element)
					if( listFixations[g][iFix] == tempCm1[g][iTemp] ){
						tempCm1[g].erase( tempCm1[g].begin() + iTemp );
						iFix ++;
					}
					//if no match, move to next element in tempCm1
					else{
						iTemp ++;
					}
				}
			
				//remove fixations from Cm2
				iFix = 0;
				iTemp = 0;
				while ( iFix < listFixations[g].size() ){
					if( listFixations[g][iFix] == tempCm2[g][iTemp] ){
						tempCm2[g].erase( tempCm2[g].begin() + iTemp );
						iFix ++;
					}
					else{
						iTemp ++;
					}
				}
			}
		}
		//*****Remove fixed mutation from all individuals and Cm
		
	}
	//*****Record and remove fixations if listFixationSize > 0
}


//*****Eddie MODDED Mar 17 2016
//Modify addMutation_V3 such that we allow for duplicate newly chosen mutations
//	These mutations will be ignored and not added to the gamete (Same as Ahmad's code)
void Popul::addMutation_V4(vector< vector<int> >& Gamete1, vector< vector<int> >& Gamete2, vector< vector<int> >& groupedLoci, int CheckForError) {

	// These variables save the average number of mutations on each chromosome
	double haploidDelMutRate, haploidBenMutRate;
    
	//haploid mutation rate = 1/2 of MUTATIONRATE (diploid mutation rate)
	haploidDelMutRate = (0.5) * MUTATIONRATE;
	haploidBenMutRate = (0.5) * U_BEN;
	
	//Number of new mutations to add to gamete 1 Cm (G1) and gamete 2 Cm (G2)
	int numNewMutG1_Del, numNewMutG1_Ben, numNewMutG2_Del, numNewMutG2_Ben;

	//Initialize random number generator, change default seed using Gamete[1].size()
	const gsl_rng_type * T;
	gsl_rng * r;
	T = gsl_rng_default;
	r = gsl_rng_alloc (T);
	gsl_rng_default_seed += time(NULL);
	gsl_rng_set(r, gsl_rng_default_seed);

	//random Poisson number of new mutations according to means calculated above
	numNewMutG1_Del = gsl_ran_poisson(r,haploidDelMutRate);
	numNewMutG1_Ben = gsl_ran_poisson(r,haploidBenMutRate);
	numNewMutG2_Del = gsl_ran_poisson(r,haploidDelMutRate);
	numNewMutG2_Ben = gsl_ran_poisson(r,haploidBenMutRate);
	
	//cout << numNewMutG1_Del << "," << numNewMutG1_Ben << "\n";

	//2D vectors to hold the chosen new mutations, each contain 3 vectors for mut group 0,1,2
	vector< vector<int> > vNewMutG1_All (3, vector<int> ());
	vector< vector<int> > vNewMutG1_Del (3, vector<int> ());
	vector< vector<int> > vNewMutG1_Ben (3, vector<int> ());
	
	vector< vector<int> > vNewMutG2_All (3, vector<int> ());
	vector< vector<int> > vNewMutG2_Del (3, vector<int> ());
	vector< vector<int> > vNewMutG2_Ben (3, vector<int> ());
		
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
	//chosenLoci record what loci has been already added to individual
	vector<int> chosenLoci;
	
	//*****Choose new mut for Gamete1
	//Choose mutations based on multinomial sampling from set of locus groups 0,1,2
	//Prevents choosing duplicate new mutations WITHIN a gamete
		
	//set countMut to 0 and clear list of chosenLoci
	//	chosenLoci are always positive integers
	//		never negative even if there are beneficials because it just tracks which loci has gained a mutation (regardless of del or ben) so I won't mutate a loci twice
	//		positive and negative integer is placed into Cm1[] and Cm2[] to indicate type of mutation
	countMut = 0;
	chosenLoci.clear();

	//choose new mutations until countMut == (numNewMutG1_Del+numNewMutG1_Ben); i.e. all ben and del mutations chosen
	while( countMut < (numNewMutG1_Del+numNewMutG1_Ben) ){
		
		//newly chosen mutation
		int tempNewMut;
		//indicate group of newly chosen mut
		int newLocus[3];
		//duplicate loci indicator, initialize as 0
		int duplicateLoci;

		//Choose new mutation, if it is a duplicate, duplicateLoci = 1, and we DO NOT add it to list of new mutations
		gsl_rng_default_seed += time(NULL);
		gsl_rng_set(r, gsl_rng_default_seed);
		//initially assume duplicate is unique
		duplicateLoci = 0;
		//reset new class of new mut
		newLocus[0] = 0; 
		newLocus[1] = 0;
		newLocus[2] = 0;
	
		//uniform random variate
		double UnifRV = gsl_ran_flat(r, 0.0, 1.0);
		
		if(countMut < numNewMutG1_Del){	
			//choose mut group 0 if fall between 0 and Bin[0]	
			if( UnifRV <= Bin[0]){
				int tempIndex = floor(gsl_ran_flat(r, -0.5, (numLoci[0]-1)+0.4999) + 0.5);
				tempNewMut = groupedLoci[0][tempIndex];
				newLocus[0] = 1;
			}
			//choose mut from group 1 if fall between Bin[0] and Bin[1]
			else if ( UnifRV > Bin[0] && UnifRV <= Bin[1]){
				int tempIndex = floor(gsl_ran_flat(r, -0.5, (numLoci[1]-1)+0.4999) + 0.5);
				tempNewMut = groupedLoci[1][tempIndex];
				newLocus[1] = 1;
			}
			//choose mut from group 2 if fall between Bin[1] and 1
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
					
		//check through list of chosenLoci to check if mutation already chosen
		vector<int>::iterator itI = chosenLoci.begin();
		while( duplicateLoci == 0 && itI < chosenLoci.end() ){
			if( tempNewMut == *itI ){
				duplicateLoci = 1;
			}
			itI++;
		}

	
		//Record mutation in vNewMutG1_Del[] if it is NOT already chosen (i.e. it doesn't exist in chosenLoci)
		//	No need to add duplicate mutation in chosenLoci, since it already exists in there
		if(duplicateLoci == 0){	
			//given new mut is non duplicate within chosenLoci, add to correct mut group; ben loci multiple by -1
			if(newLocus[0] == 1){
				if(countMut < numNewMutG1_Del){
					vNewMutG1_Del[0].push_back(tempNewMut);
				}
				else{
					vNewMutG1_Ben[0].push_back(tempNewMut * -1);
				}
			}
			else if(newLocus[1] == 1){
				if(countMut < numNewMutG1_Del){
					vNewMutG1_Del[1].push_back(tempNewMut);
				}
				else{
					vNewMutG1_Ben[1].push_back(tempNewMut * -1);
				}
			}
			else if(newLocus[2] == 1){
				if(countMut < numNewMutG1_Del){
					vNewMutG1_Del[2].push_back(tempNewMut);
				}
				else{
					vNewMutG1_Ben[2].push_back(tempNewMut * -1);
				}
			}
		
			//add tempNewMut to vector of chosenLoci_Cm1 (used to check for duplicates)
			chosenLoci.push_back(tempNewMut);
		}
		//Prints out when we sample repeat mutations, it is not added to gamete NOR resampled
		/*
		else{
			cout << "\n" << "GAM 1: CHOSEN REPEAT MUTATION OUT OF " << numNewMutG1_Del+numNewMutG1_Ben;
			if(newLocus[0] == 1){
				cout << " in loci group 0 (" << numLoci[0] << ")" << "\n";
			}
			else if(newLocus[1] == 1){
				cout << " in loci group 1(" << numLoci[1] << ")"  << "\n";
			}
			else if(newLocus[2] == 1){
				cout << " in loci group 2(" << numLoci[2] << ")"  << "\n";
			}
			
			cout << "Only added " << chosenLoci.size() << " out of " << numNewMutG1_Del << " mutations" << "\n";
		}
		*/
		
		//increase count of new mutations
		countMut++;
	}
	
	//*****Choose new mut for Gamete1
	
	
	//*****Choose new mut for Gamete2
	//set countMut to 0 and clear list of chosenLoci
	countMut = 0;
	chosenLoci.clear();

	//choose new mutations until countMut == (numNewMutG2_Del+numNewMutG2_Ben)
	while( countMut < (numNewMutG2_Del+numNewMutG2_Ben) ){
		
		//newly chosen mutation
		int tempNewMut;
		//indicate group of newly chosen mut
		int newLocus[3];
		//duplicate loci indicator, initialize as 0
		int duplicateLoci;


		gsl_rng_default_seed += time(NULL)^countMut;
		gsl_rng_set(r, gsl_rng_default_seed);
		//initially assume duplicate is unique
		duplicateLoci = 0;
		//reset new class of new mut
		newLocus[0] = 0; 
		newLocus[1] = 0;
		newLocus[2] = 0;
	
		//uniform random variate
		double UnifRV = gsl_ran_flat(r, 0.0, 1.0);
		
		if(countMut < numNewMutG2_Del){	
			if( UnifRV <= Bin[0]){
				int tempIndex = floor(gsl_ran_flat(r, -0.5, (numLoci[0]-1)+0.4999) + 0.5);
				tempNewMut = groupedLoci[0][tempIndex];
				newLocus[0] = 1;
			}
			else if ( UnifRV > Bin[0] && UnifRV <= Bin[1]){
				int tempIndex = floor(gsl_ran_flat(r, -0.5, (numLoci[1]-1)+0.4999) + 0.5);
				tempNewMut = groupedLoci[1][tempIndex];
				newLocus[1] = 1;
			}
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
		
					
		//check through list of chosenLoci to check if mutation already exists 
		vector<int>::iterator itI = chosenLoci.begin();
		while( duplicateLoci == 0 && itI < chosenLoci.end() ){
			if( tempNewMut == *itI ){
				duplicateLoci = 1;
			}
			itI++;
		}

		
				
		//Record mutation in vNewMutG1_Del[] if it is NOT already chosen (i.e. it doesn't exist in chosenLoci)
		//	No need to add duplicate mutation in chosenLoci, since it already exists in there
		if(duplicateLoci == 0){	
			if(newLocus[0] == 1){
				if(countMut < numNewMutG2_Del){
					vNewMutG2_Del[0].push_back(tempNewMut);
				}
				else{
					vNewMutG2_Ben[0].push_back(tempNewMut * -1);
				}
			}
			else if(newLocus[1] == 1){
				if(countMut < numNewMutG2_Del){
					vNewMutG2_Del[1].push_back(tempNewMut);
				}
				else{
					vNewMutG2_Ben[1].push_back(tempNewMut * -1);
				}
			}
			else if(newLocus[2] == 1){
				if(countMut < numNewMutG2_Del){
					vNewMutG2_Del[2].push_back(tempNewMut);
				}
				else{
					vNewMutG2_Ben[2].push_back(tempNewMut * -1);
				}
			}
		
			//add tempNewMut to vector of chosenLoci (used to check for duplicates)
			chosenLoci.push_back(tempNewMut);
		}
		
		//increase count of new mutations
		countMut++;
		
	}
	//*****Choose new mut for Gamete2
	
	//*****Construct and sort list of new mutations (Del are positive, ben are negative, but sorted by absolute value)
	for(int i = 0; i < 3; i++){
		//combine del and ben mutations
		vNewMutG1_All[i].insert(vNewMutG1_All[i].end(), vNewMutG1_Del[i].begin(), vNewMutG1_Del[i].end());
		vNewMutG1_All[i].insert(vNewMutG1_All[i].end(), vNewMutG1_Ben[i].begin(), vNewMutG1_Ben[i].end());
		vNewMutG2_All[i].insert(vNewMutG2_All[i].end(), vNewMutG2_Del[i].begin(), vNewMutG2_Del[i].end());
		vNewMutG2_All[i].insert(vNewMutG2_All[i].end(), vNewMutG2_Ben[i].begin(), vNewMutG2_Ben[i].end());
		//sort combined list by absolute value
		sort(vNewMutG1_All[i].begin(), vNewMutG1_All[i].end(), compareAbsolute);
		sort(vNewMutG2_All[i].begin(), vNewMutG2_All[i].end(), compareAbsolute);
	}
	//*****Construct and sort list of new mutations (Del are positive, ben are negative, but sorted by absolute value using the compareAbsolute functor)
	
	//View list of new mut chosen
	/*
	cout << "Add:\n";
	cout << numNewMutG1_Del << "," << numNewMutG2_Del << "," << numNewMutG1_Ben << "," << numNewMutG2_Ben << "\n";
	cout << "vNewMutG1_All:" << "\n";
	PrintVecInt2D(vNewMutG1_All);
	cout << "vNewMutG2_All:" << "\n";
	PrintVecInt2D(vNewMutG2_All);
	*/

	
	//*****Add new mutations to existing parental chromosomes
	//Ignores new mutation if already exist in Gamete; change signs if new mut is opposite sign of one existent in gamete
	if( CheckForError == 0 ){
		addMutToGamete(vNewMutG1_All, Gamete1);
		addMutToGamete(vNewMutG2_All, Gamete2);
	}
	//check for error in adding new mut by comparing fast algorithm (used above) to brute force method that takes more time
	else if( CheckForError == 1 ){
		addMutToGameteCheck(vNewMutG1_All, Gamete1, RUNTIME_ERRORS);
		addMutToGameteCheck(vNewMutG2_All, Gamete2, RUNTIME_ERRORS);
	}
	//*****Add new mutations to existing parental chromosomes
	
    //free memory; prevent memory leaks
	gsl_rng_free(r);
}
//*****Eddie MODDED Mar 17 2016


//Calc mean and var num of deleterious mutations
void Popul::calcDelMut(){
	//initialize with zeros
	for(int i = 0; i < 2; i++){
		populTotalDel[i] = 0;
		for(int j = 0; j < 3; j++){
			populDel[i][j] = 0;
		}
	}
	
	//store sum squares, var for mut group and for ind must be calc separately
	double sumSquaresGroup[3] = {};
	double sumSquaresIndiv = 0;
	
	//Get values from all inds
	for(int i = 0; i < POPULATION; i++){
		int temp1 = 0;
		
		for(int j = 0; j < 3; j++){
			//temp references to delMut[j] of ind i
			int& temp2 = thePop[i].delMut[j];
			
			temp1 += temp2;
			
			populDel[0][j] += temp2;
			sumSquaresGroup[j] += pow(temp2, 2.0);
		}
		sumSquaresIndiv += pow( temp1, 2.0 );
	}
	
	//Calc mean and var for each class
	for(int j = 0; j < 3; j++){
		populDel[0][j] /= POPULATION;
		populDel[1][j] = sumSquaresGroup[j]/POPULATION - pow( populDel[0][j], 2.0 );
	}
	//Calc mean and var for individual
	populTotalDel[0] = populDel[0][0] + populDel[0][1] + populDel[0][2];
	populTotalDel[1] = sumSquaresIndiv/POPULATION - pow( populTotalDel[0], 2.0 );
}


//Calc mean and var num of beneficial mutations
void Popul::calcBenMut(){
	//initialize with zeros
	for(int i = 0; i < 2; i++){
		populTotalBen[i] = 0;
		for(int j = 0; j < 3; j++){
			populBen[i][j] = 0;
		}
	}
	
	//store sum squares
	double sumSquaresGroup[3] = {};
	double sumSquaresIndiv = 0;
	//Get values from all inds
	for(int i = 0; i < POPULATION; i++){
		int temp1 = 0;
		
		for(int j = 0; j < 3; j++){
			//temp references to benMut[j] of ind i
			int& temp2 = thePop[i].benMut[j];
			temp1 += temp2;
			
			populBen[0][j] += temp2;
			sumSquaresGroup[j] += pow(temp2, 2.0);
		}
		sumSquaresIndiv += pow( temp1, 2.0 );
	}
	
	//Calc mean and var for each class
	for(int j = 0; j < 3; j++){
		populBen[0][j] /= POPULATION;
		populBen[1][j] = sumSquaresGroup[j]/POPULATION - pow( populBen[0][j], 2.0 );
	}
	//Calc mean and var for individual
	populTotalBen[0] = populBen[0][0] + populBen[0][1] + populBen[0][2];
	populTotalBen[1] = sumSquaresIndiv/POPULATION - pow( populTotalBen[0], 2.0 );
}


//Calc mean and var number of deleterious homozygous and heterozygous loci
void Popul::calcDelHeteroHomo(){

	//initialize with zeros
	for(int i = 0; i < 2; i++){
		
		populTotalDelHomozygous[i] = 0;
		populTotalDelHeterozygous[i] = 0;
	
		for(int j = 0; j < 3; j++){
			populDelHomozygous[i][j] = 0;
			populDelHeterozygous[i][j] = 0;
		}
	}

	//store sum squares
	double sumSquaresGroupHomo[3] = {};
	double sumSquaresIndivHomo = 0;	
	double sumSquaresGroupHetero[3] = {};
	double sumSquaresIndivHetero = 0;
	
	
	//Get values from all inds
	for(int i = 0; i < POPULATION; i++){
		int tempHomo1 = 0;
		int tempHetero1 = 0;
		
		for(int j = 0; j < 3; j++){
			//temp references to delHomo[j], delHetero[i] of ind i
			int& tempHomo2 = thePop[i].delHomo[j];
			int& tempHetero2 = thePop[i].delHetero[j];
			
			tempHomo1 += tempHomo2;
			tempHetero1 += tempHetero2;
			
			populDelHomozygous[0][j] += tempHomo2;
			sumSquaresGroupHomo[j] += pow(tempHomo2, 2.0);
			
			populDelHeterozygous[0][j] += tempHetero2;
			sumSquaresGroupHetero[j] += pow(tempHetero2, 2.0);
		}
		sumSquaresIndivHomo += pow( tempHomo1, 2.0 );
		sumSquaresIndivHetero += pow( tempHetero1, 2.0 );
	}
	
	//Calc mean and var for each class
	for(int j = 0; j < 3; j++){
		populDelHomozygous[0][j] /= POPULATION;
		populDelHomozygous[1][j] = sumSquaresGroupHomo[j]/POPULATION - pow( populDelHomozygous[0][j], 2.0 );
		
		populDelHeterozygous[0][j] /= POPULATION;
		populDelHeterozygous[1][j] = sumSquaresGroupHetero[j]/POPULATION - pow( populDelHeterozygous[0][j], 2.0 );
	}
	
	//Calc mean and var for individual
	populTotalDelHomozygous[0] = populDelHomozygous[0][0] + populDelHomozygous[0][1] + populDelHomozygous[0][2];
	populTotalDelHomozygous[1] = sumSquaresIndivHomo/POPULATION - pow( populTotalDelHomozygous[0], 2.0 );
	
	populTotalDelHeterozygous[0] = populDelHeterozygous[0][0] + populDelHeterozygous[0][1] + populDelHeterozygous[0][2];
	populTotalDelHeterozygous[1] = sumSquaresIndivHetero/POPULATION - pow( populTotalDelHeterozygous[0], 2.0 );
}


//Error checking function to make certain addMutToGamete creates the correct vector of mutations
//	Compares algorithm in addMutToGamete to one that is slower but always add new mutations to correct spot in gamete
//	Increase numError if resulting mut vectors from the two algorithms differ at any position
void Popul::addMutToGameteCheck(vector< vector<int> >& vMut, vector< vector<int> >& vGamete, int& numError){
	
	//vGameteCheck will add vMut to vGamete in slow algorithm, should be exact same as fast algorithm
	vector< vector<int> > vGameteCheck = vGamete;
	
	//fast algorithm (replicate of addMutToGamete)
	for(int i = 0; i < 3; i++){
		if (vGamete[i].size() == 0){
			vGamete[i] = vMut[i];
		}
		else{
			int iG = 0;
			int iM = 0;
				
			while( iM < vMut[i].size() ){
				if( abs(vMut[i][iM]) < abs(vGamete[i][iG]) ){
					vGamete[i].insert( vGamete[i].begin() + iG, vMut[i][iM] );			
					iM++;		
				}
				else if( abs(vMut[i][iM]) > abs(vGamete[i][iG]) ){
					if( iG == (vGamete[i].size() - 1) ){
						vGamete[i].insert( vGamete[i].end(), vMut[i].begin() + iM, vMut[i].end() );
						iM = vMut[i].size();
					}			
				}
				else{	
					if ( (vMut[i][iM] * vGamete[i][iG]) < 0 ){
						vGamete[i][iG] *= -1;
					}
					iM++;			
				}
				if(iG < (vGamete[i].size() - 1)){
					iG++;
				}
			}
		}
	
	}
	
	
	//slow algorithm, search from beginning of gamete vector each time a mut needs to be added
	for(int i = 0; i < 3; i++){
		if (vGameteCheck[i].size() == 0){
			vGameteCheck[i] = vMut[i];
		}
		//else vGameteCheck non-empty, need to add new mut and keep ascending order
		else{
			//index for mutation vector set to 0, will keep track of focal element of vectors 
			int iM = 0;
			
			while( iM < vMut[i].size() ){
		
				for(int iG = 0; iG < vGameteCheck[i].size(); iG++){
				
					if( abs(vMut[i][iM]) < abs(vGameteCheck[i][iG]) ){
						vGameteCheck[i].insert( vGameteCheck[i].begin() + iG, vMut[i][iM] );
						break;	
					}
					else if( abs(vMut[i][iM]) > abs(vGameteCheck[i][iG]) ){
						if( iG == (vGameteCheck[i].size() - 1) ){
							vGameteCheck[i].insert( vGameteCheck[i].end(), vMut[i].begin() + iM, vMut[i].end() );
							iM = vMut.size();
							break;
						}	
					}
					else{	
						if ( (vMut[i][iM] * vGameteCheck[i][iG]) < 0 ){
							vGameteCheck[i][iG] *= -1;
						}
						break;
					}
				}
				iM++;
			}
		}
	}
	
	//*****Nov18
	//sort vGameteCheck in ascending order by absolute value (just in case algorithm above is incorrect in inserting mutations under ascending order)
	for(int i = 0; i < 3; i++){
		sort(vGameteCheck[i].begin(), vGameteCheck[i].end(), compareAbsolute);
	}
	//*****Nov18
	
			
	//loop through the 3 vectors within vGamete and vGameteCheck
	for(int i = 0; i < 3; i++){
		//Increase error count if size of gamete vector in the two algorithms are unequal
		if( vGamete[i].size() != vGameteCheck[i].size() ){
			numError++;
			cout << vGamete[i].size() << " != " << vGameteCheck[i].size() << "\n";
		}
		//Increase error count if the values of the two gamete vectors are not in identical order
		else{
			for(int j = 0; j < vGamete[i].size(); j++){
				if(vGamete[i][j] != vGameteCheck[i][j]){
					numError++;
				}
			}
		}
	}
}

//*****Eddie MODDED Mar 7 2016
//Check that the loci in vecFixations are actually fixed by counting if ALL Cm of ALL individuals have mutations at those loci
//Input: reference to vecFixations 
void Popul::checkFixations(vector< vector<int> >& vecFixations){
	
	//Count of number of matches to loci in vecFixations among ALL Cm of ALL individuals
	//fixationMatches should equal vecFixationsize*2*POPULATION after finish counting
	fixationMatches = 0;
	
	//Record number of loci in vecFixations
	vecFixationSize = vecFixations[0].size() + vecFixations[1].size() + vecFixations[2].size();
	
	//loop through all n individuals					
	for (int n = 0; n < POPULATION; n++){
		
		//listTemp REFERENCES to individual n's Cm1 for comparing to listFixation
		vector< vector<int> >& listTemp1 = thePop[n].Cm1;
		vector< vector<int> >& listTemp2 = thePop[n].Cm2;
		
		//Cm1: loop through all three mutation groups
		for(int g = 0; g < 3; g++){

			//index element in for vecFixations[g][]
			int iFix = 0;

			//loop through list of mutated loci in group g
			for(int h = 0; h < listTemp1[g].size(); h++){
								
				//if iFix at end of vecFixations[g], break out of for loop (looked through all of vecFixations[g] already)
				if( vecFixations[g].size() == iFix){
					break;
				}
				if( vecFixations[g][iFix] == listTemp1[g][h] ){
					fixationMatches++;
					iFix++;
				}
			}
		}
	
		//Cm2: loop through all three mutation groups
		for(int g = 0; g < 3; g++){
		
			//index element in for vecFixations[g][]
			int iFix = 0;

			//loop through list of mutated loci in group g
			for(int h = 0; h < listTemp2[g].size(); h++){
				
				//if iFix at end of vecFixations[g], break out of for loop (looked through all of vecFixations[g] already)
				if( vecFixations[g].size() == iFix){
					break;
				}
				if( vecFixations[g][iFix] == listTemp2[g][h] ){
					fixationMatches++;
					iFix++;
				}
			}
		}	
	}
	//cout << "check: " << vecFixationSize << "," << fixationMatches << "\n";
}
//*****Eddie MODDED Mar 7 2016














//*****Eddie MODDED Aug 2016
//NEW constructors and function for running NEUTRAL model (no fitness affecting loci or selection)
//	for runModel2

//Constructor: When an object of type "Popul" is created, it automatically creates its first generation
Popul::Popul() {

	//Construct individuals and save in thePop
	for(int i=0; i<POPULATION; i++) {
		gsl_rng_default_seed += time(NULL)^i;
		
		//create indiv newIndiv
		Indiv newIndiv;
		
		//save newIndiv in thePop (vector of indiv)
		thePop.push_back(newIndiv);
	}
	
	//Update population measurements
	calcRecombRate();
	
	//*****Eddie MODDED Aug, 30, 2016
	//initialize neutralGameteIndex with POPULATION,3,2,2 elements in dimension 1,2,3,4, respectively
	neutralGameteIndex.resize(POPULATION);
	for(int n = 0; n < POPULATION; n++){
		neutralGameteIndex[n].resize(3);
		
		for(int l = 0; l < 3; l++){
			neutralGameteIndex[n][l].resize(2);
			
			for(int a = 0; a < 2; a++){
				neutralGameteIndex[n][l][a].resize(2);
				
				//initialize the parental indices as -1; will give error unless the value is changed to 0 or 1
				neutralGameteIndex[n][l][a][0]=-1;
				neutralGameteIndex[n][l][a][1]=-1;
			}
		}
	}
	//*****Eddie MODDED Aug, 30, 2016
		
}

//Create new generation of individuals by calling makeOffspring
//thePop replaced by individuals in thePopOffspring and population measurements are updated
//Requires groupedLoci since involves adding new mutations (requires groupedLoci)
void Popul::makeNewGeneration_New(){
	
	//*****Eddie MODDED Aug, 30, 2016
	//re-initialize the parental indices as -1; will give error unless the value is changed to 0 or 1
	for(int n = 0; n < POPULATION; n++){		
		for(int l = 0; l < 3; l++){			
			for(int a = 0; a < 2; a++){				
				neutralGameteIndex[n][l][a][0]=-1;
				neutralGameteIndex[n][l][a][1]=-1;
			}
		}
	}
	//*****Eddie MODDED Aug, 30, 2016
	
	//Initialize summary statistic of parent ID to 0 each generation	
	ID_P1 = 0;
	ID_P2 = 0;
	
	//create POPULATION number of offspring stored into thePopOffspring
	for(int i = 0; i < POPULATION; i++) {
		makeOffspring_New(i);
	}

	//Calc average parental ID
	ID_P1/=POPULATION;
	ID_P2/=POPULATION;
	
	//*****thePop now updated as thePopOffspring, creating the next generation
	thePop = thePopOffspring;
	
	//Clear thePopOffspring, extremely important when making offspring for the generation after this
	thePopOffspring.clear();
	
	//Update population measurements
	calcRecombRate();
}



//select index of first parent for producing offspring
//return index of parent that pass selection
int Popul::pickFirstParent_New(){
	
	//index of indiv in thePop that passes selection
	int indexParent;

	//Initialize RNG	
	const gsl_rng_type * T;
	gsl_rng * r;
	T = gsl_rng_taus;
	r = gsl_rng_alloc (T);
	gsl_rng_default_seed += time(NULL);
	gsl_rng_set (r, gsl_rng_default_seed);
	
	//random pick index from 0 to POPULATION-1 (for population with size N, it is indexed from 0 to N-1)
	indexParent = floor(gsl_ran_flat(r, -0.5, (POPULATION-1)+0.4999) + 0.5);
	
	//free memory; prevent memory leaks
	gsl_rng_free(r);
	
	//return index of indiv that passed fitness check
	return indexParent;
}



//select index of parent for outcrossing, utilizes pollen discounted fitness (indivPDFitness)
//return index of parent that pass selection based on pollen discounted fitness
int Popul::pickOutcrossParent_New( int indexParent1 ){
	
	//index of indiv in thePop that passes selection
	int indexParent2;
	
	//Initialize RNG	
	const gsl_rng_type * T;
	gsl_rng * r;
	T = gsl_rng_taus;
	r = gsl_rng_alloc (T);
	gsl_rng_default_seed += time(NULL);
	gsl_rng_set (r, gsl_rng_default_seed);

		
	//random pick index from 0 to POPULATION-1 that is != indexParent1
	do{
		indexParent2 = floor(gsl_ran_flat(r, -0.5, (POPULATION-1)+0.4999) + 0.5);
	}while(indexParent2 == indexParent1);
	
	//free memory; prevent memory leaks
	gsl_rng_free(r);
	
	//return index of indiv that passed pollen discounted fitness check
	return indexParent2;
}


//Create a SINGLE offspring individual
//1. Select parents based on fitness and polen discounted fitness if outcrossing
//2. Parents create gametes with recombination rate of populRecombRate[0]
//3. Add new mutations to gametes
//4. Mutate modifier loci if required
//5. Create Indiv Offspring given gamete1, gamete2 and modifier loci values
//6. Store Offspring in thePopOffspring vector
void Popul::makeOffspring_New(int n){
	
	//Initialize RNG	
	const gsl_rng_type * T;
	gsl_rng * r;
	T = gsl_rng_taus;
	r = gsl_rng_alloc (T);
	gsl_rng_default_seed += time(NULL);
	gsl_rng_set (r, gsl_rng_default_seed);
	
	//index of parent 1 and parent 2 within thePop
	int indexP1, indexP2;
	
	//*****Selection: select parents based on their individual fitness relative to the population highest fitness
	//P1 outcrosses if selfing rate < Uniform random variable [0,1], P2 chosen based on pollen discount fitness
	//P1 selfs if selfing rate >= Uniform random variable [0,1], P2 = P1; (i.e. always selfes if selfing rate == 1)
	
	//select parent 1 based on indivFitness
	indexP1 = pickFirstParent_New();

	//decide if P1 will self or outcross
	double UnifRV = gsl_ran_flat(r, 0, 1);
	
	//Initialize OffspringID
	int OffspringID = 0;
	
	//P1 outcross with P2 != P1 (P2 selected using pollen discounted fitness)
	if( UnifRV > thePop[indexP1].indivSelfingRate ) {
		indexP2 = pickOutcrossParent_New(indexP1);
		OffspringID = 1;
	}
	//P1 self, so P2 = P1
	else {
		indexP2 = indexP1;
		OffspringID = 0;
	}
	
	//Add up parent1 and parent2 IDs
	ID_P1 += thePop[indexP1].ID;
	ID_P2 += thePop[indexP2].ID;
	//*****Selection: select parents based on their individual fitness relative to population highest fitness
	
	
	//*****Create gamete for parent 1,2; individuals are thePop[indexP1] and thePop[indexP2]
	double modifierRecomb_gamete1;
	double modifierSelfing_gamete1;
	double modifierRecomb_gamete2;
	double modifierSelfing_gamete2;
	
	//indicate which Cm of P1,P2 to inherit neutral allele
	int neutral_P1Allele[3] = {};
	int neutral_P2Allele[3] = {};
		
	//makeGamete for parent 1,2 and store info into appropriate variables
	thePop[indexP1].makeGamete_New( modifierRecomb_gamete1, modifierSelfing_gamete1, neutral_P1Allele );
	thePop[indexP2].makeGamete_New( modifierRecomb_gamete2, modifierSelfing_gamete2, neutral_P2Allele );
	//*****Create gamete for parent 1,2
	
	
	//*****Create offspring using constructor of Indiv
	Indiv Offspring( modifierRecomb_gamete1, modifierRecomb_gamete2, modifierSelfing_gamete1, modifierSelfing_gamete2);
	//*****Create offspring using constructor of Indiv	
	
	Offspring.ID = OffspringID;
	
	//*****Store Offspring in thePopOffspring
	//Assumes thePopOffspring has been cleared in the previous generation
	thePopOffspring.push_back( Offspring );
	//*****Store Offspring in thePopOffspring
	
	//free memory; prevent memory leaks
	gsl_rng_free(r);
	
	
	//*****Record neutral indices
	//For offspring n's allele 0 for all three neutral loci
	neutralGameteIndex[n][0][0][0]=indexP1;
	neutralGameteIndex[n][0][0][1]=neutral_P1Allele[0];
	neutralGameteIndex[n][1][0][0]=indexP1;
	neutralGameteIndex[n][1][0][1]=neutral_P1Allele[1];
	neutralGameteIndex[n][2][0][0]=indexP1;
	neutralGameteIndex[n][2][0][1]=neutral_P1Allele[2];
	
	//For offspring n's allele 1 for all three neutral loci
	neutralGameteIndex[n][0][1][0]=indexP2;
	neutralGameteIndex[n][0][1][1]=neutral_P2Allele[0];
	neutralGameteIndex[n][1][1][0]=indexP2;
	neutralGameteIndex[n][1][1][1]=neutral_P2Allele[1];
	neutralGameteIndex[n][2][1][0]=indexP2;
	neutralGameteIndex[n][2][1][1]=neutral_P2Allele[2];

	/*Display neutralGameteIndex info
	cout << "Print neutralGameteIndex:\n";
	cout << neutralGameteIndex[n][0][0][0] << "," << neutralGameteIndex[n][0][0][1] << "\n";
	cout << neutralGameteIndex[n][1][0][0] << "," << neutralGameteIndex[n][1][0][1] << "\n";
	cout << neutralGameteIndex[n][2][0][0] << "," << neutralGameteIndex[n][2][0][1] << "\n";
	cout << neutralGameteIndex[n][0][1][0] << "," << neutralGameteIndex[n][0][1][1] << "\n";
	cout << neutralGameteIndex[n][1][1][0] << "," << neutralGameteIndex[n][1][1][1] << "\n";
	cout << neutralGameteIndex[n][2][1][0] << "," << neutralGameteIndex[n][2][1][1] << "\n";
	*/
	//*****Record neutral indices	
}

//*****Eddie MODDED Aug 2016










//*****Eddie MODDED Apr 20 2017
//Modified addMutation_V4 to ignore adding mutations to fluctuating or constant loci

//Same as addMutation_V4 except ignore adding mutations to constant loci (loci in groupedLoci[2])
//	still allows mutations to occur but mutations not added to Cm1, Cm2
void Popul::addMutation_V4_IGNORE1(vector< vector<int> >& Gamete1, vector< vector<int> >& Gamete2, vector< vector<int> >& groupedLoci, int CheckForError) {

	// These variables save the average number of mutations on each chromosome
	double haploidDelMutRate, haploidBenMutRate;
    
	//haploid mutation rate = 1/2 of MUTATIONRATE (diploid mutation rate)
	haploidDelMutRate = (0.5) * MUTATIONRATE;
	haploidBenMutRate = (0.5) * U_BEN;

	//Number of new mutations to add to gamete 1 Cm (G1) and gamete 2 Cm (G2)
	int numNewMutG1_Del, numNewMutG1_Ben, numNewMutG2_Del, numNewMutG2_Ben;

	//Initialize random number generator, change default seed using Gamete[1].size()
	const gsl_rng_type * T;
	gsl_rng * r;
	T = gsl_rng_default;
	r = gsl_rng_alloc (T);
	gsl_rng_default_seed += time(NULL);
	gsl_rng_set(r, gsl_rng_default_seed);

	//random Poisson number of new mutations according to means calculated above
	numNewMutG1_Del = gsl_ran_poisson(r,haploidDelMutRate);
	numNewMutG1_Ben = gsl_ran_poisson(r,haploidBenMutRate);
	numNewMutG2_Del = gsl_ran_poisson(r,haploidDelMutRate);
	numNewMutG2_Ben = gsl_ran_poisson(r,haploidBenMutRate);
	

	//2D vectors to hold the chosen new mutations, each contain 3 vectors for mut group 0,1,2
	vector< vector<int> > vNewMutG1_All (3, vector<int> ());
	vector< vector<int> > vNewMutG1_Del (3, vector<int> ());
	vector< vector<int> > vNewMutG1_Ben (3, vector<int> ());
	
	vector< vector<int> > vNewMutG2_All (3, vector<int> ());
	vector< vector<int> > vNewMutG2_Del (3, vector<int> ());
	vector< vector<int> > vNewMutG2_Ben (3, vector<int> ());
		
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
	//chosenLoci record what loci has been already added to individual
	vector<int> chosenLoci;
	
	//*****Choose new mut for Gamete1
	//Choose mutations based on multinomial sampling from set of locus groups 0,1,2
	//Prevents choosing duplicate new mutations WITHIN a gamete
		
	//set countMut to 0 and clear list of chosenLoci
	countMut = 0;
	chosenLoci.clear();

	//choose new mutations until countMut == (numNewMutG1_Del+numNewMutG1_Ben); i.e. all ben and del mutations chosen
	while( countMut < (numNewMutG1_Del+numNewMutG1_Ben) ){
		
		//newly chosen mutation
		int tempNewMut;
		//indicate group of newly chosen mut
		int newLocus[3];
		//duplicate loci indicator, initialize as 0
		int duplicateLoci;

		//Choose new mutation, if it is a duplicate, duplicateLoci = 1, and we DO NOT add it to list of new mutations
		gsl_rng_default_seed += time(NULL);
		gsl_rng_set(r, gsl_rng_default_seed);
		//initially assume duplicate is unique
		duplicateLoci = 0;
		//reset new class of new mut
		newLocus[0] = 0; 
		newLocus[1] = 0;
		newLocus[2] = 0;
	
		//uniform random variate
		double UnifRV = gsl_ran_flat(r, 0.0, 1.0);
	
		//choose mut group 0 if fall between 0 and Bin[0]	
		if( UnifRV <= Bin[0]){
			int tempIndex = floor(gsl_ran_flat(r, -0.5, (numLoci[0]-1)+0.4999) + 0.5);
			tempNewMut = groupedLoci[0][tempIndex];
			newLocus[0] = 1;
		}
		//choose mut from group 1 if fall between Bin[0] and Bin[1]
		else if ( UnifRV > Bin[0] && UnifRV <= Bin[1]){
			int tempIndex = floor(gsl_ran_flat(r, -0.5, (numLoci[1]-1)+0.4999) + 0.5);
			tempNewMut = groupedLoci[1][tempIndex];
			newLocus[1] = 1;
		}

		
					
		//check through list of chosenLoci to check if mutation already chosen
		vector<int>::iterator itI = chosenLoci.begin();
		while( duplicateLoci == 0 && itI < chosenLoci.end() ){
			if( tempNewMut == *itI ){
				duplicateLoci = 1;
			}
			itI++;
		}

	
		//Record mutation in vNewMutG1_Del[] if it is NOT already chosen (i.e. it doesn't exist in chosenLoci)
		//	No need to add duplicate mutation in chosenLoci, since it already exists in there
		if(duplicateLoci == 0){	
			//given new mut is non duplicate within chosenLoci, add to correct mut group; ben loci multiple by -1
			if(newLocus[0] == 1){
				if(countMut < numNewMutG1_Del){
					vNewMutG1_Del[0].push_back(tempNewMut);
				}
				else{
					vNewMutG1_Ben[0].push_back(tempNewMut * -1);
				}
			}
			else if(newLocus[1] == 1){
				if(countMut < numNewMutG1_Del){
					vNewMutG1_Del[1].push_back(tempNewMut);
				}
				else{
					vNewMutG1_Ben[1].push_back(tempNewMut * -1);
				}
			}

		
			//add tempNewMut to vector of chosenLoci_Cm1 (used to check for duplicates)
			chosenLoci.push_back(tempNewMut);
		}
		
		//increase count of new mutations
		//	increase count even if constant loci was chosen (we still treat it as a mutation event)
		countMut++;
	}
	
	//*****Choose new mut for Gamete1
	
	
	//*****Choose new mut for Gamete2
	//set countMut to 0 and clear list of chosenLoci
	countMut = 0;
	chosenLoci.clear();

	//choose new mutations until countMut == (numNewMutG2_Del+numNewMutG2_Ben)
	while( countMut < (numNewMutG2_Del+numNewMutG2_Ben) ){
		
		//newly chosen mutation
		int tempNewMut;
		//indicate group of newly chosen mut
		int newLocus[3];
		//duplicate loci indicator, initialize as 0
		int duplicateLoci;


		gsl_rng_default_seed += time(NULL)^countMut;
		gsl_rng_set(r, gsl_rng_default_seed);
		//initially assume duplicate is unique
		duplicateLoci = 0;
		//reset new class of new mut
		newLocus[0] = 0; 
		newLocus[1] = 0;
		newLocus[2] = 0;
	
		//uniform random variate
		double UnifRV = gsl_ran_flat(r, 0.0, 1.0);
	
		//choose mut group 0 if fall between 0 and Bin[0]	
		if( UnifRV <= Bin[0]){
			int tempIndex = floor(gsl_ran_flat(r, -0.5, (numLoci[0]-1)+0.4999) + 0.5);
			tempNewMut = groupedLoci[0][tempIndex];
			newLocus[0] = 1;
		}
		//choose mut from group 1 if fall between Bin[0] and Bin[1]
		else if ( UnifRV > Bin[0] && UnifRV <= Bin[1]){
			int tempIndex = floor(gsl_ran_flat(r, -0.5, (numLoci[1]-1)+0.4999) + 0.5);
			tempNewMut = groupedLoci[1][tempIndex];
			newLocus[1] = 1;
		}
		
					
		//check through list of chosenLoci to check if mutation already exists 
		vector<int>::iterator itI = chosenLoci.begin();
		while( duplicateLoci == 0 && itI < chosenLoci.end() ){
			if( tempNewMut == *itI ){
				duplicateLoci = 1;
			}
			itI++;
		}

		
				
		//Record mutation in vNewMutG1_Del[] if it is NOT already chosen (i.e. it doesn't exist in chosenLoci)
		//	No need to add duplicate mutation in chosenLoci, since it already exists in there
		if(duplicateLoci == 0){	
			if(newLocus[0] == 1){
				if(countMut < numNewMutG2_Del){
					vNewMutG2_Del[0].push_back(tempNewMut);
				}
				else{
					vNewMutG2_Ben[0].push_back(tempNewMut * -1);
				}
			}
			else if(newLocus[1] == 1){
				if(countMut < numNewMutG2_Del){
					vNewMutG2_Del[1].push_back(tempNewMut);
				}
				else{
					vNewMutG2_Ben[1].push_back(tempNewMut * -1);
				}
			}
		
			//add tempNewMut to vector of chosenLoci (used to check for duplicates)
			chosenLoci.push_back(tempNewMut);
		}
		
		//increase count of new mutations
		//	increase count even if constant loci was chosen (we still treat it as a mutation event)
		countMut++;
		
	}
	//*****Choose new mut for Gamete2
	
	//*****Construct and sort list of new mutations (Del are positive, ben are negative, but sorted by absolute value)
	for(int i = 0; i < 3; i++){
		//combine del and ben mutations
		vNewMutG1_All[i].insert(vNewMutG1_All[i].end(), vNewMutG1_Del[i].begin(), vNewMutG1_Del[i].end());
		vNewMutG1_All[i].insert(vNewMutG1_All[i].end(), vNewMutG1_Ben[i].begin(), vNewMutG1_Ben[i].end());
		vNewMutG2_All[i].insert(vNewMutG2_All[i].end(), vNewMutG2_Del[i].begin(), vNewMutG2_Del[i].end());
		vNewMutG2_All[i].insert(vNewMutG2_All[i].end(), vNewMutG2_Ben[i].begin(), vNewMutG2_Ben[i].end());
		//sort combined list by absolute value
		sort(vNewMutG1_All[i].begin(), vNewMutG1_All[i].end(), compareAbsolute);
		sort(vNewMutG2_All[i].begin(), vNewMutG2_All[i].end(), compareAbsolute);
	}
	//*****Construct and sort list of new mutations (Del are positive, ben are negative, but sorted by absolute value using the compareAbsolute functor)
	
	
	//*****Add new mutations to existing parental chromosomes
	//Ignores new mutation if already exist in Gamete; change signs if new mut is opposite sign of one existent in gamete
	if( CheckForError == 0 ){
		addMutToGamete(vNewMutG1_All, Gamete1);
		addMutToGamete(vNewMutG2_All, Gamete2);
	}
	//check for error in adding new mut by comparing fast algorithm (used above) to brute force method that takes more time
	else if( CheckForError == 1 ){
		addMutToGameteCheck(vNewMutG1_All, Gamete1, RUNTIME_ERRORS);
		addMutToGameteCheck(vNewMutG2_All, Gamete2, RUNTIME_ERRORS);
	}
	//*****Add new mutations to existing parental chromosomes

	
    //free memory; prevent memory leaks
	gsl_rng_free(r);
}











//Same as addMutation_V4 except ignore adding mutations to fluctuating loci (loci in groupedLoci[0], groupedLoci[1])
//	still allows mutations to occur but mutations not added to Cm1, Cm2
void Popul::addMutation_V4_IGNORE2(vector< vector<int> >& Gamete1, vector< vector<int> >& Gamete2, vector< vector<int> >& groupedLoci, int CheckForError) {

	// These variables save the average number of mutations on each chromosome
	double haploidDelMutRate, haploidBenMutRate;
    
	//haploid mutation rate = 1/2 of MUTATIONRATE (diploid mutation rate)
	haploidDelMutRate = (0.5) * MUTATIONRATE;
	haploidBenMutRate = (0.5) * U_BEN;

	//Number of new mutations to add to gamete 1 Cm (G1) and gamete 2 Cm (G2)
	int numNewMutG1_Del, numNewMutG1_Ben, numNewMutG2_Del, numNewMutG2_Ben;

	//Initialize random number generator, change default seed using Gamete[1].size()
	const gsl_rng_type * T;
	gsl_rng * r;
	T = gsl_rng_default;
	r = gsl_rng_alloc (T);
	gsl_rng_default_seed += time(NULL);
	gsl_rng_set(r, gsl_rng_default_seed);

	//random Poisson number of new mutations according to means calculated above
	numNewMutG1_Del = gsl_ran_poisson(r,haploidDelMutRate);
	numNewMutG1_Ben = gsl_ran_poisson(r,haploidBenMutRate);
	numNewMutG2_Del = gsl_ran_poisson(r,haploidDelMutRate);
	numNewMutG2_Ben = gsl_ran_poisson(r,haploidBenMutRate);
	

	//2D vectors to hold the chosen new mutations, each contain 3 vectors for mut group 0,1,2
	vector< vector<int> > vNewMutG1_All (3, vector<int> ());
	vector< vector<int> > vNewMutG1_Del (3, vector<int> ());
	vector< vector<int> > vNewMutG1_Ben (3, vector<int> ());
	
	vector< vector<int> > vNewMutG2_All (3, vector<int> ());
	vector< vector<int> > vNewMutG2_Del (3, vector<int> ());
	vector< vector<int> > vNewMutG2_Ben (3, vector<int> ());
		
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
	//chosenLoci record what loci has been already added to individual
	vector<int> chosenLoci;
	
	//*****Choose new mut for Gamete1
	//Choose mutations based on multinomial sampling from set of locus groups 0,1,2
	//Prevents choosing duplicate new mutations WITHIN a gamete
		
	//set countMut to 0 and clear list of chosenLoci
	countMut = 0;
	chosenLoci.clear();

	//choose new mutations until countMut == (numNewMutG1_Del+numNewMutG1_Ben); i.e. all ben and del mutations chosen
	while( countMut < (numNewMutG1_Del+numNewMutG1_Ben) ){
		
		//newly chosen mutation
		int tempNewMut;
		//indicate group of newly chosen mut
		int newLocus[3];
		//duplicate loci indicator, initialize as 0
		int duplicateLoci;

		//Choose new mutation, if it is a duplicate, duplicateLoci = 1, and we DO NOT add it to list of new mutations
		gsl_rng_default_seed += time(NULL);
		gsl_rng_set(r, gsl_rng_default_seed);
		//initially assume duplicate is unique
		duplicateLoci = 0;
		//reset new class of new mut
		newLocus[0] = 0; 
		newLocus[1] = 0;
		newLocus[2] = 0;
	
		//uniform random variate
		double UnifRV = gsl_ran_flat(r, 0.0, 1.0);
	
		//choose mut from group 2 if UnifRV fall in Bin[2] (i.e. > Bin[1])
		if( UnifRV > Bin[1]){
			int tempIndex = floor(gsl_ran_flat(r, -0.5, (numLoci[2]-1)+0.4999) + 0.5);
			tempNewMut = groupedLoci[2][tempIndex];
			newLocus[2] = 1;
		}
		
					
		//check through list of chosenLoci to check if mutation already chosen
		vector<int>::iterator itI = chosenLoci.begin();
		while( duplicateLoci == 0 && itI < chosenLoci.end() ){
			if( tempNewMut == *itI ){
				duplicateLoci = 1;
			}
			itI++;
		}

	
		//Record mutation in vNewMutG1_Del[] if it is NOT already chosen (i.e. it doesn't exist in chosenLoci)
		//	No need to add duplicate mutation in chosenLoci, since it already exists in there
		if(duplicateLoci == 0){	
			//given new mut is non duplicate within chosenLoci, add to correct mut group; ben loci multiple by -1
			if(newLocus[2] == 1){
				if(countMut < numNewMutG1_Del){
					vNewMutG1_Del[2].push_back(tempNewMut);
				}
				else{
					vNewMutG1_Ben[2].push_back(tempNewMut * -1);
				}
			}
		
			//add tempNewMut to vector of chosenLoci_Cm1 (used to check for duplicates)
			chosenLoci.push_back(tempNewMut);
		}
		
		//increase count of new mutations
		//	increase count even if constant loci was chosen (we still treat it as a mutation event)
		countMut++;
	}
	
	//*****Choose new mut for Gamete1
	
	
	//*****Choose new mut for Gamete2
	//set countMut to 0 and clear list of chosenLoci
	countMut = 0;
	chosenLoci.clear();

	//choose new mutations until countMut == (numNewMutG2_Del+numNewMutG2_Ben)
	while( countMut < (numNewMutG2_Del+numNewMutG2_Ben) ){
		
		//newly chosen mutation
		int tempNewMut;
		//indicate group of newly chosen mut
		int newLocus[3];
		//duplicate loci indicator, initialize as 0
		int duplicateLoci;


		gsl_rng_default_seed += time(NULL)^countMut;
		gsl_rng_set(r, gsl_rng_default_seed);
		//initially assume duplicate is unique
		duplicateLoci = 0;
		//reset new class of new mut
		newLocus[0] = 0; 
		newLocus[1] = 0;
		newLocus[2] = 0;
	
		//uniform random variate
		double UnifRV = gsl_ran_flat(r, 0.0, 1.0);
	
		//choose mut from group 2 if UnifRV fall in Bin[2] (i.e. > Bin[1])
		if( UnifRV > Bin[1]){
			int tempIndex = floor(gsl_ran_flat(r, -0.5, (numLoci[2]-1)+0.4999) + 0.5);
			tempNewMut = groupedLoci[2][tempIndex];
			newLocus[2] = 1;
		}
		
					
		//check through list of chosenLoci to check if mutation already exists 
		vector<int>::iterator itI = chosenLoci.begin();
		while( duplicateLoci == 0 && itI < chosenLoci.end() ){
			if( tempNewMut == *itI ){
				duplicateLoci = 1;
			}
			itI++;
		}

		
				
		//Record mutation in vNewMutG1_Del[] if it is NOT already chosen (i.e. it doesn't exist in chosenLoci)
		//	No need to add duplicate mutation in chosenLoci, since it already exists in there
		if(duplicateLoci == 0){	
			if(newLocus[2] == 1){
				if(countMut < numNewMutG2_Del){
					vNewMutG2_Del[2].push_back(tempNewMut);
				}
				else{
					vNewMutG2_Ben[2].push_back(tempNewMut * -1);
				}
			}
		
			//add tempNewMut to vector of chosenLoci (used to check for duplicates)
			chosenLoci.push_back(tempNewMut);
		}
		
		//increase count of new mutations
		//	increase count even if constant loci was chosen (we still treat it as a mutation event)
		countMut++;
		
	}
	//*****Choose new mut for Gamete2
	
	//*****Construct and sort list of new mutations (Del are positive, ben are negative, but sorted by absolute value)
	for(int i = 0; i < 3; i++){
		//combine del and ben mutations
		vNewMutG1_All[i].insert(vNewMutG1_All[i].end(), vNewMutG1_Del[i].begin(), vNewMutG1_Del[i].end());
		vNewMutG1_All[i].insert(vNewMutG1_All[i].end(), vNewMutG1_Ben[i].begin(), vNewMutG1_Ben[i].end());
		vNewMutG2_All[i].insert(vNewMutG2_All[i].end(), vNewMutG2_Del[i].begin(), vNewMutG2_Del[i].end());
		vNewMutG2_All[i].insert(vNewMutG2_All[i].end(), vNewMutG2_Ben[i].begin(), vNewMutG2_Ben[i].end());
		//sort combined list by absolute value
		sort(vNewMutG1_All[i].begin(), vNewMutG1_All[i].end(), compareAbsolute);
		sort(vNewMutG2_All[i].begin(), vNewMutG2_All[i].end(), compareAbsolute);
	}
	//*****Construct and sort list of new mutations (Del are positive, ben are negative, but sorted by absolute value using the compareAbsolute functor)
	
	//*****Add new mutations to existing parental chromosomes
	//Ignores new mutation if already exist in Gamete; change signs if new mut is opposite sign of one existent in gamete
	if( CheckForError == 0 ){
		addMutToGamete(vNewMutG1_All, Gamete1);
		addMutToGamete(vNewMutG2_All, Gamete2);
	}
	//check for error in adding new mut by comparing fast algorithm (used above) to brute force method that takes more time
	else if( CheckForError == 1 ){
		addMutToGameteCheck(vNewMutG1_All, Gamete1, RUNTIME_ERRORS);
		addMutToGameteCheck(vNewMutG2_All, Gamete2, RUNTIME_ERRORS);
	}
	//*****Add new mutations to existing parental chromosomes

	
    //free memory; prevent memory leaks
	gsl_rng_free(r);
}
//*****Eddie MODDED Apr 20 2017



