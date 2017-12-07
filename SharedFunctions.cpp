#include "SharedFunctions.h"

//Given 2D int vector, will print out all contents of vVec[i][] line by line by i
void PrintVecInt2D(vector< vector<int> >& vVec)
{
    for(int i = 0; i < vVec.size(); i++){
		for(int j = 0; j < vVec[i].size(); j++){
			cout << vVec[i][j] << " ";
		}
		cout << endl;
    }
}

void PrintVecInt(vector<int>& vVec)
{
    for(int i = 0; i < vVec.size(); i++){
			cout << vVec[i] << " ";
    }
    cout << endl;
}


// a helper function to compare two integers by their absolute values
bool compareAbsolute(int x, int y) {
	return abs(x)<abs(y);
}


//Calc mean and var individual given a vector of integers
void calcMeanVar_Double(vector<double>& vecDouble, double &vecMean, double &vecVar){
	//Initialize mean and var
	vecMean = 0;
	vecVar = 0;
	
	double sumSquares = 0;
	int vecSize = vecDouble.size();
	
	//calc sum and sum squares of vecInteger
	for(int i = 0; i < vecSize; i++){
		vecMean += vecDouble[i];
		sumSquares += pow( vecDouble[i], 2.0 );
	}
	
	//calc mean and var
	vecMean /= vecSize;
	vecVar = sumSquares/vecSize - pow( vecMean, 2.0 );
}

//calc and return correlation between vec1 and vec2
double calcCorrelation_Double(vector<double>& vec1, vector<double>& vec2){
	//size of vec1 and vec2
	int size1 = vec1.size();
	int size2 = vec2.size();
	
	double cor12=0;
	//calc correlation of vectors or equal size
	if(size1==size2){
		double mean1, mean2, var1, var2;
		double cov12;
		
		//calc mean and var of vec1 and vec2
		calcMeanVar_Double(vec1,mean1,var1);
		calcMeanVar_Double(vec2,mean2,var2);
		//calc covariance
		for(int i=0; i < size1; i++){
			cov12 += (vec1[i]-mean1)*(vec2[i]-mean2);
		}
		cov12/=size1;
		//calc correlation 
		cor12= cov12/(sqrt(var1)*sqrt(var2));
		//return correlation
		return(cor12);
		
	}
	else{
		//return 0 is vec1 and vec2 different size
		cout << "ERROR: vector 1 and vector 2 have different sizes." << "\n";
		return (0);
	}
}