#include "InputPar.h"


int InitialErrorCheck (string SaveDir, int NumSim, int NumSimNames){

	int ErrorCount = 0;

	//file name starts as "SaveDir"
	string ErrorFileName(SaveDir);
	ErrorFileName.append("InitialErrorFile");
	ErrorFileName.append(".txt");
	
	ofstream ErrorFile;
	ErrorFile.open (ErrorFileName.c_str(), ios::out);
		
		//Check SaveDir exists
		if(SaveDir.size() == 0){
			ErrorFile << "ERROR: Save directory NOT specified" << "\n";
			ErrorCount++;
		}
		else{
			ErrorFile << "No error, save directory specified" << "\n";
		}
	
		//Check NumSim and NumSimNames match
		if(NumSim != NumSimNames){
			ErrorFile << "ERROR: NumSim != NumSimNames." << "\n";
			ErrorCount++;
		}
		else{
			ErrorFile << "No error, equal number of simulations and number of simulations names" << "\n";
		}
		
		if(ErrorCount > 0){
			ErrorFile << "Simulation stopped. Did NOT proceed to check parameters." << "\n";
		}
		
	ErrorFile.close();

	//Return count of errors
	return ErrorCount;
	
}



//Check for a number of different errors in parameter value combinations
//Returns number of errors, simulation will not start if there is one or more errors
int ParameterErrorCheck(string SaveDir, string FileName){

	int ErrorCount = 0;

	//file name starts as "SaveDir"
	string ErrorFileName(SaveDir);
	ErrorFileName.append(FileName);
	ErrorFileName.append("_ParErrorFile");
	ErrorFileName.append(".txt");
	
	ofstream ErrorFile;
	ErrorFile.open (ErrorFileName.c_str(), ios::out | ios::app);
	
		ErrorFile << "Parameter Error Check File for: " << FileName << "\n\n";
		
		//Check proportion/frequency parameters within range of 0 and 1
		if(INITIALSELFING > 1 || F_AC > 1 || PHI > 1 || PSI > 1 || T0 > 1 || T1 > 1 || (S*DELTA_S0_A) > 1 || (S*DELTA_S0_B) > 1 || (S*DELTA_S1_A) > 1 || (S*DELTA_S1_B) > 1 || (S*DELTA_S2_C) > 1){
			ErrorFile << "ERROR: Proportion/frequency parameters not within 0 and 1" << "\n";
			ErrorCount++;
		}
		else{
			ErrorFile << "No error in proportion/frequency parameters" << "\n";
		}

		//Check indicator variables are binary
		if(STOCHASTIC_FLUCTUATION < 0 || STOCHASTIC_FLUCTUATION > 2){
			ErrorFile << "ERROR: indicator variables not binary" << "\n";
			ErrorCount++;
		}
		else{
			ErrorFile << "No error for indicator variables" << "\n";
		}

		//Check deterministic fluctuation values
		if( (T0*TAU)+T1 > 1 || MAXNUMOFGENERATIONS < PERIOD){
			ErrorCount++;
			ErrorFile << "ERROR in deterministic fluctuations: (EVEN_A*TAU)+ODD_A > 1 or MAXNUMOFGENERATIONS < PERIOD" << "\n";
		}
		else{
			ErrorFile << "No error for generating deterministic fluctuations" << "\n";
		}

		//Check stochastic fluctuation values
		//Assign probabilities of each environment occurring
		double EnvProb[4];
		EnvProb[0] = (PHI)*(PSI) + EPSILON;
		EnvProb[1] = (PHI)*(1-PSI) - EPSILON;
		EnvProb[2] = (1-PHI)*(PSI) - EPSILON;
		EnvProb[3] = (1-PHI)*(1-PSI) + EPSILON;

		if( EnvProb[0] > 1 || EnvProb[1] > 1 || EnvProb[2] > 1 || EnvProb[3] > 1 || EnvProb[0] < 0 || EnvProb[1] < 0 || EnvProb[2] < 0 || EnvProb[3] < 0){
			ErrorCount++;
			ErrorFile << "ERROR stochastic fluctuations: Probabilities of environment occurrence NOT between 0 and 1" << "\n";
		}
		else{
			ErrorFile << "No error for generating stochastic fluctuations" << "\n";
		}

		//Check if file name is default
		if(MODELFILENAME == "DEFAULT_MODEL_FILE_NAME"){
			ErrorFile << "ERROR: Using default file name" << "\n";
			ErrorCount++;
		}
		else{
			ErrorFile << "No error in file name, using user defined file name" << "\n";
		}
		

		//Output number of errors/warnings
		ErrorFile << "\n" << "Number of errors: " << ErrorCount << "\n";
	
		if(ErrorCount == 0){
			ErrorFile << "No errors, yay!" << "\n";
		}
		else if(ErrorCount == 8){
			ErrorFile << "Achievement Unlocked: Maximum number of errors!";
		}
		
		
		ErrorFile 	<< 	"\n";
		ErrorFile	<<	"MODELFILENAME " << MODELFILENAME << "\n";
		ErrorFile	<<	"NUMOFREPS " << NUMOFREPS << "\n";
		ErrorFile	<<	"NUMOFGENERATIONS " << NUMOFGENERATIONS << "\n";
		ErrorFile	<<	"MAXNUMOFGENERATIONS " << MAXNUMOFGENERATIONS << "\n";
		ErrorFile	<<	"L " << L << "\n";
		ErrorFile	<<	"POPULATION " << POPULATION << "\n";
		ErrorFile	<<	"MUTATIONRATE " << MUTATIONRATE << "\n";
		ErrorFile	<<	"H " << H << "\n";
		ErrorFile	<<	"S " << S << "\n";
		ErrorFile	<<	"INITIALSELFING " << INITIALSELFING << "\n";
		ErrorFile	<<	"INITIALRECOMBINATION " << INITIALRECOMBINATION << "\n";
		ErrorFile	<<	"STOCHASTIC_FLUCTUATION " << STOCHASTIC_FLUCTUATION << "\n";
		ErrorFile	<<	"PERIOD " << PERIOD << "\n";
		ErrorFile	<<	"T0 " << T0 << "\n";
		ErrorFile	<<	"T1 " << T1 << "\n";
		ErrorFile	<<	"TAU " << TAU << "\n";
		ErrorFile	<<	"F_AC " << F_AC << "\n";
		ErrorFile	<<	"PHI " << PHI << "\n";
		ErrorFile	<<	"PSI " << PSI << "\n";
		ErrorFile	<<	"EPSILON " << EPSILON << "\n";
		ErrorFile	<<	"ETA " << ETA << "\n";
		ErrorFile	<<	"SAVE_DATA_INTERVAL " << SAVE_DATA_INTERVAL << "\n";
		ErrorFile	<<	"K " << K << "\n";
		ErrorFile	<<	"SELFINGMUTATIONRATE " << SELFINGMUTATIONRATE << "\n";
		ErrorFile	<<	"SIGMAS " << SIGMAS << "\n";
		ErrorFile	<<	"FINALSELFINGMUTATIONRATE " << FINALSELFINGMUTATIONRATE << "\n";
		ErrorFile	<<	"FINALSIGMAS " << FINALSIGMAS << "\n";
		ErrorFile	<<	"RECOMBMUTATIONRATE " << RECOMBMUTATIONRATE << "\n";
		ErrorFile	<<	"SIGMAR " << SIGMAR << "\n";
		ErrorFile	<<	"DELTA_S0_A " << DELTA_S0_A << "\n";
		ErrorFile	<<	"DELTA_S0_B " << DELTA_S0_B << "\n";
		ErrorFile	<<	"DELTA_S1_A " << DELTA_S1_A << "\n";
		ErrorFile	<<	"DELTA_S1_B " << DELTA_S1_B << "\n";
		ErrorFile	<<	"DELTA_S2_C " << DELTA_S2_C << "\n";
		ErrorFile	<<	"P_CONSTANT " << P_CONSTANT << "\n";
		ErrorFile	<<	"P_FLUC_G0 " << P_FLUC_G0 << "\n";
	
	ErrorFile.close();
	
	//Return count of errors
	return ErrorCount;
}




//Generate Example Parameter File 
//ExampleParameterFile.txt will NOT allow simulation to run because it does NOT specify a saving directory (-dir has no following string value)
void WriteExampleParameterFile(){
	
	//Make an Example parameter file
	string FileName = "ExampleParameterFile.txt";
	ofstream ExampleFile;
	//ios::out simply sets the file as an output file
	//add .c_str() to make string FILENAME because .open only takes const char* as parameter
	ExampleFile.open (FileName.c_str(), ios::out);
	ExampleFile << "     \n";
	ExampleFile << "-ng 2500\n";
	ExampleFile << "-mng 2500\n";
	ExampleFile << "     -l 5000#Can have spaces on the left of the line, since it will be trimmed, best avoid this\n";	
	ExampleFile << "-n 10000   #Can have spaces on the right of the line, since it will be trimmed, best avoid this\n";
	ExampleFile << "-u 1.0\n";
	ExampleFile << "-h 0.25\n";
	ExampleFile << "-s 0.05\n";
	ExampleFile << "-is 1\n";
	ExampleFile << "-ir 5\n";
	ExampleFile << "-sto 1\n";
	ExampleFile << "-pe 80\n";
	ExampleFile << "-t0 0.5\n";
	ExampleFile << "-t1 0.5\n";
	ExampleFile << "-tau 0\n";
	ExampleFile << "-f 0.95\n";
	ExampleFile << "-phi 0.5\n";
	ExampleFile << "-psi 0.5\n";
	ExampleFile << "-eps 0.25\n";
	ExampleFile << "-eta 0\n";
	ExampleFile << "-xdat 0\n";
	ExampleFile << "-k 0\n";
	ExampleFile << "-ius 0\n";
	ExampleFile << "-iss 0\n";
	ExampleFile << "-fss 0\n";
	ExampleFile << "-fus 0\n";
	ExampleFile << "-iur 0\n";
	ExampleFile << "-isr 0\n";
	ExampleFile << "-mfn ExampleName1 ExampleName2 ExampleName3 ExampleName4\n";
	ExampleFile << "-dir\n";
	ExampleFile << "-mv 1\n";
	ExampleFile << "-d0a 1";
	ExampleFile << "-d0b 0";
	ExampleFile << "-d1a 1";
	ExampleFile << "-d1b 0";
	ExampleFile << "-d2c 0.5";
	ExampleFile << "-prc 0.4";
	ExampleFile << "-prf0 0.5";
	ExampleFile << "#No blank line at last line";
	ExampleFile.close();
	
}

	
//Function to read keys, values and filenames from parameter text file with directory/name ParFileName
//Will pass reference to vector<double> Keys, vector<vector<double>> Values, vector<string> SimNames to store parameters
//Will pass reference to # of simulations to run (NumSims) (i.e. different par values, NOT replicates) and # of simulation names (NumSimNames)
//NumSims must be EQUAL to NumSimNames, otherwise there are unequal number of simulations and simulation names
//NumSims and NumSimNames will be passed in to ErrorCheck() function later to make sure they are equal

//sKeys will NOT contain "-mfn" since this is key to filenames (thus nValues will also NOT contain filenames)
//Filenames of simulations will be stored in SimNames instead
void ReadParameterFile (vector<string> &sKeys, vector< vector<double> > &nValues, vector<string> &SimNames, string &SaveDir, int &NumSims, int &NumSimNames, string ParFileName){

	int 	NumSimsCount = 1;
	int 	NumSimNamesCount = 0;
	
	//default Input file directory/name; its the file we just created above	
	string INPUTFILENAME = ParFileName;

	ifstream InputFile;
	InputFile.open(INPUTFILENAME.c_str(), ios::in);
	while (!InputFile.eof()){
        
        std::string strInput;
		vector<string> RemoveComment;
		vector<string> splitStr;

		//get current line and split into substrings based on delimiter ' ' (i.e. white space)
        getline(InputFile, strInput);
        
        //remove comments from  line
        //only need to keep RemoveComment[0] (i.e. code BEFORE comment)
        split(strInput, '#', RemoveComment);
        
		//left and right trim white spaces of RemoveComment[0] (RemoveComment[0] are contents of strInput without comments)
        const char* ws = " \t\n\r\f\v";
        rtrim(RemoveComment[0],ws);
        ltrim(RemoveComment[0],ws);

        //skip current iteration if RemoveComment[0] is empty after removing comment and white spaces
        if(RemoveComment[0].size() < 1){
        //cout << "Skip blank line" << "\n";
        continue;
        }
        //if line not empty, there is useful code, split by white space
        else{
		split(RemoveComment[0], ' ', splitStr);
		}
		
		////if splitStr[0] == "-dir", it means splitStr[1] is the directory to save result files
		//Save directory stored into SaveDirectory (there should only be ONE save directory)
		if(splitStr[0] == "-dir"){
			
			//store value of -dir if splitStr.size > 1 (i.e. splitStr[1] exists and has some value)
			if(	splitStr.size() > 1){		
			SaveDir = splitStr[1];
			}
		}
		//if splitStr[0] == "-mfn", it means splitStr[1+] give model file names
		//Combine the model file name with SaveDirectory and store in SimNames
		else if(splitStr[0] == "-mfn"){                            
			//for loop starts at element 1 of splitStr (skip splitStr[0], which is "-mfn")
			for(int i = 1; i < splitStr.size(); i++){                                                                			
				SimNames.push_back( splitStr[i] );
				
				//count number of simulation names
				NumSimNamesCount++;		
			}
		}
		//if splitStr[0] != "-mfn", splitStr contains parameters for simulation
		else{
			//This section records the Keys and Values after removing comments and spaces
			//Keys and Values are stored as different elements of splitStr
			//Keys are splitStr[0], any element past 0 are Values for that Key
			
			//keys are element 0 of splitStr
			sKeys.push_back(splitStr[0]);
		
			//convert subset of splitStr (element 1 to end) into type double
			//tempNum will hold converted value of each element in splitStr
			//tempVec will hold all tempNum (i.e. tempVec = splitStr converted to double)
			vector<double> tempVec;
			double tempNum;
			vector<string>::iterator itI;                                         
			//for loop starts at element 1 (splitStr.begin()+1)
			for(itI = splitStr.begin() + 1; itI < splitStr.end(); itI++){                                                                
				ConvertStrToDouble(*itI,tempNum);
			
				tempVec.push_back(tempNum);
													  
			}
			//save tempVec in to nValues
			nValues.push_back(tempVec);
				
			//number of par combo that will be made
			NumSimsCount*=tempVec.size();
		}

    }	
	
	InputFile.close();

	//Pass by reference the count of # sims and # sim names
	//If parameter file was formatted correctly NumSims = NumSimNames
	NumSims = NumSimsCount;
	NumSimNames = NumSimNamesCount;

}




//Given MainVec is a 2-D vector with MainVec.size() rows and MainVec[i].size() columns
//Creates all possible combinations of elements by column going down each row
//records a completed column to 2-D vector Save (passed by reference) 
//should have Cumulative-Product of MainVec[i].size() number of vectors at the end
void MakeParameterList(vector< vector<double> >& MainVec, size_t index, vector<double> SubVec, vector< vector<double> >& Save){
	
	//Called at last recursion after going through all rows
	if(index >= MainVec.size()){
		
		//print SubVec here to show the list at the final recursion
		//PrintVecDouble(&SubVec);
		
		//record vector Sum in to Save; Save is 2-D vector
		Save.push_back(SubVec);
		
		return;
	}
	//When not at last row, add element in current row and continue to next row
	//when current row has more than one element, then need to create additional recursions (i.e. for loop need to loop more than once)
	else{
	
		vector<double>::iterator itI;
		for(itI=MainVec[index].begin(); itI < MainVec[index].end(); itI++)                   //itI initialized with vVec->begin()
		{                                                                   					//loops as long as itI is smaller than vVec->end()			
			//add element in current row to SubVec
			SubVec.push_back(*itI);
			MakeParameterList(MainVec,index+1,SubVec,Save);
			
			if(MainVec[index].size()>1){
				//pop_back to remove last element
				//allow last element to be REPLACED by next element in MainVec[index]
				//if statement is redundant
				//	since if MainVec[index].size() ==1, loop exits even in pop_back incorrectly occurs                                         
				SubVec.pop_back();
			}
		}
	}
}


//Assign parameters in ParList[SimNumber][] to the simulation global parameter variables
void AssignParameterValues ( vector<string>& vKeys, vector< vector<double> >& vParList, int SimNumber, string SaveDir, string FileName ){
	
	for(int j=0; j < vKeys.size(); j++){
		
		//cout << vKeys[j] << " " << vParList[SimNumber][j] << "\n";

		if( vKeys[j] == "-nr" ){
			NUMOFREPS = vParList[SimNumber][j];
		}
		else if( vKeys[j] == "-ng" ){
			NUMOFGENERATIONS = vParList[SimNumber][j];
		}
		else if( vKeys[j] == "-mng" ){
			MAXNUMOFGENERATIONS = vParList[SimNumber][j];
		}
		else if( vKeys[j] == "-l" ){
			L = vParList[SimNumber][j];
		}
		else if( vKeys[j] == "-n" ){
			POPULATION = vParList[SimNumber][j];
		}
		else if( vKeys[j] == "-u" ){
			MUTATIONRATE = vParList[SimNumber][j];
		}
		else if( vKeys[j] == "-h" ){
			H = vParList[SimNumber][j];
		}
		else if( vKeys[j] == "-s" ){
			S = vParList[SimNumber][j];
		}
		else if( vKeys[j] == "-k" ){
			K = vParList[SimNumber][j];
		}
		else if( vKeys[j] == "-is" ){
			INITIALSELFING = vParList[SimNumber][j];
		}
		else if( vKeys[j] == "-ius" ){
			SELFINGMUTATIONRATE = vParList[SimNumber][j];
		}
		else if( vKeys[j] == "-iss" ){
			SIGMAS = vParList[SimNumber][j];
		}
		else if( vKeys[j] == "-fus" ){
			FINALSELFINGMUTATIONRATE = vParList[SimNumber][j];
		}
		else if( vKeys[j] == "-fss" ){
			FINALSIGMAS = vParList[SimNumber][j];
		}
		else if( vKeys[j] == "-ir" ){
			INITIALRECOMBINATION = vParList[SimNumber][j];
		}
		else if( vKeys[j] == "-iur" ){
			RECOMBMUTATIONRATE = vParList[SimNumber][j];
		}
		else if( vKeys[j] == "-isr" ){
			SIGMAR = vParList[SimNumber][j];
		}
		else if( vKeys[j] == "-sto" ){
			STOCHASTIC_FLUCTUATION = vParList[SimNumber][j];
		}
		else if( vKeys[j] == "-pe" ){
			PERIOD = vParList[SimNumber][j];
		}
		else if( vKeys[j] == "-t0" ){
			T0 = vParList[SimNumber][j];
		}
		else if( vKeys[j] == "-t1" ){
			T1 = vParList[SimNumber][j];
		}
		else if( vKeys[j] == "-tau" ){
			TAU = vParList[SimNumber][j];
		}
		else if( vKeys[j] == "-f" ){
			F_AC = vParList[SimNumber][j];
		}
		else if( vKeys[j] == "-phi" ){
			PHI = vParList[SimNumber][j];
		}
		else if( vKeys[j] == "-psi" ){
			PSI = vParList[SimNumber][j];
		}
		else if( vKeys[j] == "-eps" ){
			EPSILON = vParList[SimNumber][j];
		}
		else if( vKeys[j] == "-xdat" ){
			SAVE_DATA_INTERVAL = vParList[SimNumber][j];
		}
		else if( vKeys[j] == "-d0a" ){
			DELTA_S0_A = vParList[SimNumber][j];
		}
		else if( vKeys[j] == "-d0b" ){
			DELTA_S0_B = vParList[SimNumber][j];
		}
		else if( vKeys[j] == "-d1a" ){
			DELTA_S1_A = vParList[SimNumber][j];
		}
		else if( vKeys[j] == "-d1b" ){
			DELTA_S1_B = vParList[SimNumber][j];
		}
		else if( vKeys[j] == "-d2c" ){
			DELTA_S2_C = vParList[SimNumber][j];
		}
		else if( vKeys[j] == "-prc" ){
			P_CONSTANT = vParList[SimNumber][j];
		}
		else if( vKeys[j] == "-prf0" ){
			P_FLUC_G0 = vParList[SimNumber][j];
		}
		else if( vKeys[j] == "-chfx" ){
			CHECK_FIXATION = vParList[SimNumber][j];
		}
		else if( vKeys[j] == "-xmin" ){
			SAVE_DATA_MIN_GEN = vParList[SimNumber][j];
		}
		else if( vKeys[j] == "-qneu" ){
			Q_NEUTRAL_LOCI = vParList[SimNumber][j];
		}
		else if( vKeys[j] == "-eta" ){
			ETA = vParList[SimNumber][j];
		}
		else if( vKeys[j] == "-sl" ){
			SELF_LENGTH = vParList[SimNumber][j];
		}
		else if( vKeys[j] == "-igl" ){
			IGNORE_LOCI = vParList[SimNumber][j];
		}
		else if( vKeys[j] == "-pfix" ){
			PFIX_MODEL = vParList[SimNumber][j];
		}
		else if( vKeys[j] == "-pfixh" ){
			PFIX_H = vParList[SimNumber][j];
		}
		else if( vKeys[j] == "-pfixs" ){
			PFIX_S = vParList[SimNumber][j];
		}
		else if( vKeys[j] == "-pfixk" ){
			PFIX_K = vParList[SimNumber][j];
		}
		else if( vKeys[j] == "-ub" ){
			U_BEN = vParList[SimNumber][j];
		}
		else if( vKeys[j] == "-sb" ){
			S_BEN = vParList[SimNumber][j];
		}
		else if( vKeys[j] == "-hb" ){
			H_BEN = vParList[SimNumber][j];
		}
		//Any unrecognized key will be report in ErrorFile
		else{
		
			//file name starts as "SaveDir"
			string ErrorFileName(SaveDir);
			ErrorFileName.append(FileName);
			ErrorFileName.append("_ParErrorFile");
			ErrorFileName.append(".txt");
			
			ofstream ErrorFile;
			ErrorFile.open (ErrorFileName.c_str(), ios::out | ios::app);
			
			ErrorFile << "WARNING: Unrecognized key: " << vKeys[j] << "\n";
			
			ErrorFile.close();
		}
		
		//Determine whether simulation is set to have fluctuating environments
		//These are the subset of conditions that I use to create const-sel, there are other possibilities that I NEVER use
		if( (DELTA_S0_A==DELTA_S0_B && DELTA_S1_A==DELTA_S1_B) || (STOCHASTIC_FLUCTUATION == 1 && PHI==1 && PSI==1) || (STOCHASTIC_FLUCTUATION == 0 && T0==1 && T1==1) ){
			ENV_FLUC_INDICATOR = 0;
		}
		else{
			ENV_FLUC_INDICATOR = 1;
		}
	}
}




//Prints Keys and parameter list OutputFileName, adds Round to end of filename
//output as Comma Separated file
void PrintParameterListToFile(vector<string>& vKeys, vector< vector<double> >& vValues, string OutputFileName, int Round)
{
	stringstream sstr;
	sstr << Round;
	//file name starts as "OutputFileName", then sstr (=Round) is appended to it
	string FileName(OutputFileName);
	FileName.append(sstr.str());
	FileName.append(".txt");

	//ios::out simply sets the file as an output file
	//add .c_str() to make string FileName because .open only takes const char* as parameter
	ofstream OutputFile;
	OutputFile.open (FileName.c_str(), ios::out);
	
	//Output parameter Keys as top line
	for(int i = 0; i < vKeys.size(); i++){
		OutputFile << vKeys[i];
		
		//if at last element of vVec[i], put "\n", otherwise put comma
		if( i == (vKeys.size() - 1) ){
			OutputFile << "\n";	
		}
		else{
			OutputFile << ",";
		}
	}
	
	//Output parameter Values
	for(int i = 0; i < vValues.size(); i++){
		for(int j = 0; j < vValues[i].size(); j++){
			OutputFile << vValues[i][j];
			
			//if at last element of vVec[i], put "\n", otherwise put comma
			if( j == (vValues[i].size() - 1) ){
				OutputFile << "\n";	
			}
			else{
				OutputFile << ",";
			}
		}
    }

	OutputFile.close();
}	
	
	
	

//Given 2D double vector, will print out all contents of vVec[i][] line by line by i
void PrintVecDouble2D(vector< vector<double> >& vVec)
{
    for(int i = 0; i < vVec.size(); i++){
		for(int j = 0; j < vVec[i].size(); j++){
			cout << vVec[i][j] << " ";
		}
		cout << endl;
    }
}




//trims const char* t from the END of string "s"
//rtrim = right trim
void rtrim(string& s, const char* t){
    s.erase( s.find_last_not_of(t) + 1 );
}

//trims const char* t from the BEGINNING of string "s"
//ltrim = left trim
void ltrim(string& s,  const char* t){
   s.erase( 0, s.find_first_not_of(t) );
}

//input full string "s" to split/separate it at anywhere that character "c" occurs
//store the separated components of "s" in to different elements of "v" in same order as "s"
//for white spaces, typically use: const char* ws = " \t\n\r\f\v"
void split(const string& s, char c, vector<string>& v) {
	string::size_type i = 0;
	string::size_type j = s.find(c);

	if(j == string::npos){
		//cout << "No match to '" << c << "' in string. String NOT split." << "\n";
   		v.push_back(s.substr(i, s.length()));
	}
	else{
		while (j != string::npos) {
			v.push_back(s.substr(i, j-i));
			i = ++j;
			j = s.find(c, j);

			if (j == string::npos){
				v.push_back(s.substr(i, s.length()));
			}
		}
	}
}


//Attempts to convert some string Text in to double Number; gives error and Number=0 if undoable
//Number is pass by reference
void ConvertStrToDouble(string Text, double& Number){

	istringstream convert (Text);
	
	//Text is converted in Number within if statement
	//If undoable, gives error and Number = 0
	if( !(convert >> Number)){
		cout << "Error: cannot convert string to double";
		Number = 0;
	}
	convert.clear();	
}



void ConvertStrToInt(string Text, int& Number){

	istringstream convert (Text);
	
	//Text is converted in Number within if statement
	//If undoable, gives error and Number = 0
	if( !(convert >> Number)){
		cout << "Error: cannot convert string to double";
		Number = 0;
	}
	convert.clear();	
}







