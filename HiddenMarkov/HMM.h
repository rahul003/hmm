#include <iostream>
#include <string>
#include "config.h"

using namespace std;

// Type for the number of observation sequece
typedef unsigned long type_T;
// Type for the number of states
typedef unsigned short type_N;
// Type for the number of possible observations OR size of codebook
typedef unsigned int type_M;
// Type used for probabilities. Can be changed to 128 bit data type
typedef long double type_prob;

class HMM{

	// TODO : Implement taking these param in constructor ( currently defined in config.h)
	/*
	type_T T;		// Number of data in observation sequece
	type_N N;		// Number of states
	type_M M;		// Range of observation data OR size of codebook
	*/
public:
	
	
	// All array are 1 based indexed

	// The Hidden Markov Model Parameters
	// The Aux Parameters
	
	long double HMM::predict(vector <unsigned int> obs, int T);
	void calculateAlpha(int T);
	void calculateBeta(int T);
	void calculateGamma(int T);
	void calculateXi(int T);
	void viterbi(int T);
	void updateModel(int T);
	void saveModel(string fileName);
	
	template<class In> 
	void printAlpha(In & file, int T);
	template<class In> 
	void printBeta(In & file, int T);
	template<class In> 
	void printGamma(In & file, int T);
	template<class In> 
	void printXi(In & file, int T);
	template<class In> 
	void printPsi(In & file, int T);
	template<class In> 
	void printDelta(In & file, int T);
	template<class In> 
	void printA(In & file);
	template<class In> 
	void printB(In & file);
	template<class In> 
	void printPI(In & file);
	template<class In> 
	void printALL(In & file, int T);

	HMM();
	void init(int T);
	void loadModel(string filepath);
	//HMM(type_prob Ain[N+1][M+1], type_prob Bin[N+1][M+1], type_prob piin[N+1]);
	HMM(char * fileName);
	void setObservationSeq(vector <unsigned int> arr, int T);
 	void train(vector <unsigned int> observations, int T);
};