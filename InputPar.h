//Need to include "included-files.h" because it links a bunch of functions and it links "theConstants.h"
//Just to to include "included-files.h" and "theConstants.h" when constructing InputPar.0 in the MakeFile
#include "IncludeFiles.h"

int InitialErrorCheck (string SaveDir, int NumSim, int NumSimNames);

int ParameterErrorCheck (string SaveDir, string FileName);

void WriteExampleParameterFile ();

void ReadParameterFile (vector<string> &sKeys, vector< vector<double> > &nValues, vector<string> &SimNames, string &SaveDir, int &NumSims, int &NumSimNames, string ParFileName);

void MakeParameterList (vector< vector<double> >& MainVec, size_t index, vector<double> SubVec, vector< vector<double> >& Save);

void AssignParameterValues (vector<string>& vKeys, vector< vector<double> >& vParList, int SimNumber, string SaveDir, string FileName);

void PrintParameterListToFile (vector<string>& vKeys, vector< vector<double> >& vValues, string OutputFileName, int Round);


void PrintVecDouble2D (vector< vector<double> >& vVec);

void rtrim(string& s, const char* t);

void ltrim(string& s,  const char* t);

void split(const string& s, char c, vector<string>& v);

void ConvertStrToDouble(string Text, double& Number);
