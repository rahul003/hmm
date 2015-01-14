#include <iostream>
#include <cstring>
#include <iterator>
#include <fstream>
#include <vector>
#include "HMM.h"
#include <sstream>
using namespace std;
int N = 5;
int M = 32;
vector < vector <long double> > alpha;	// alpha[t][i] is the probability of observing a partial sequence of 
							// observables o1,…ot such that at time t, state qt=i
vector < vector <long double> > beta;	// beta[t][i] – the probability of observing a sequence of
							// observables ot+1,…,oT given state qt =i at time t, and model
vector < vector <long double> > gamma;
vector <vector < vector <type_prob> > > xi;
vector <type_N> qStar;
vector <type_M> O;// Observations
	// For Viterbi Algorithm
vector < vector <unsigned short> > psi;
vector <vector <type_prob> > delta;    // delta[t][i] – the probability of the most probable path ending in state qt=i
vector <vector <type_prob> > A;
vector <vector <type_prob> > B;
vector <type_prob> pi;
// The state transition probability matrix, aij = P(qt+1 = j|qt = i)
// Observation probability distribution, bj(k) = P(ot = k|qt = j)  i ≤ k ≤ M
// The initial state distribution
type_prob prob_observ_given_model;
type_prob pStar;
HMM* componentmodels[10];
HMM* finalmodels[50];

string toString(int number)
{
   stringstream ss;//create a stringstream
   ss << number;//add number to the stream
   return ss.str();//return a string with the contents of the stream
}

HMM::HMM(){
	//T = num;
	//cout<<num;
	prob_observ_given_model = 0.0;
	pStar = 0.0;
	for(int i=0;i<N;i++)
	{
		vector <type_prob> temp;
		for(int j=0;j<M;j++)
		{
			temp.push_back(0);
		}
		A.push_back(temp);
	}
	for(int i=0;i<N;i++)
	{
		vector <type_prob> temp;
		for(int j=0;j<M;j++)
		{
			temp.push_back(0);
		}
		B.push_back(temp);
	}
	for(int i=0;i<N;i++)
	{
		pi.push_back(0);
	}

	//A.resize(N+1);
	//for(int i=0;i<M+1;i++)
	//{A[i].resize(M+1);}
	//B.resize(N+1);
	//for(int i=0;i<M+1;i++)
	//{B[i].resize(M+1);}
	//pi.resize(N+1);
	
	type_prob oneByN = 1.0 / N;
	type_prob oneByM = 1.0 / M;
	for(int i=0 ; i<N ; i++){
            if (i==0)
				pi[i] = 1;
			else
				pi[i]=0;
            for (int j = 0; j <N ; ++j)
            {
            	if (i==j && i!=N-1)
					A[i][j] =0.8;
				else if (i==j-1)
					A[i][j]=0.2;
				else if (i==N-1 && j==N-1)
					A[i][j]=1;
				else
					A[i][j]=0;
            }
            for (int j = 0; j <M ; ++j)
            {
            	B[i][j] = oneByM;
            }
    }
    
}

void HMM::init(int T)
{
	alpha.clear();
	beta.clear();
	gamma.clear();
	delta.clear();
	psi.clear();
	xi.clear();
	qStar.clear();
	for(int i=0;i<T;i++)
	{
		vector <type_prob> temp;
		for(int j=0;j<N;j++)
		{
			temp.push_back(0);
		}
		alpha.push_back(temp);
	}

	for(int i=0;i<T;i++)
	{
		vector <type_prob> temp;
		for(int j=0;j<N;j++)
		{
			temp.push_back(0);
		}
		beta.push_back(temp);
	}

	for(int i=0;i<T;i++)
	{
		vector <type_prob> temp;
		for(int j=0;j<N;j++)
		{
			temp.push_back(0);
		}
		gamma.push_back(temp);
	}

	for(int i=0;i<T;i++)
	{
		vector <type_prob> temp;
		for(int j=0;j<N;j++)
		{
			temp.push_back(0);
		}
		delta.push_back(temp);
	}

	
	for(int i=0;i<T;i++)
	{
		vector <type_N> temp;
		for(int j=0;j<N;j++)
		{
			temp.push_back(0);
		}
		psi.push_back(temp);
	}

	//vector<vector<vector<type_prob>>> vec1(T+1, vector <vector<type_prob> >(N+1, vector<type_prob>(N+1)));
    for(int i = 0; i < T; i++)
	  {
		vector < vector < type_prob > > w;
		xi.push_back( w );
		for(int j = 0; j < N; j++)
		{
		  vector <type_prob> v;
		  xi[i].push_back( v );
		  for(int k = 0; k < N; k++)
		  {
			xi[i][j].push_back(0);
		  }
		}
	  }
	
	for (int i = 0; i < T; ++i)
    {
        qStar.push_back(0);
    }

	for(int i=0;i<N;i++)
	{
		for(int j=0;j<M;j++)
		{
			if(B[i][j]==0)
				B[i][j]=1E-20;
		}
	}

}
/*
HMM::HMM(type_prob Ain[N+1][M+1], type_prob Bin[N+1][M+1], type_prob piin[N+1])
{
	prob_observ_given_model = 0.0;
	pStar = 0.0;
	
	type_prob oneByN = 1.0 / N;
	type_prob oneByM = 1.0 / M;
	for(int i=0 ; i<N ; i++){
		pi[i] = piin[i];
            for (int j = 0; j <N ; ++j)
            {
            	A[i][j] = Ain[i][j];
            }
            for (int j = 0; j <M ; ++j)
            {
            	B[i][j] = Bin[i][j];
            }
    }
    
}*/
void HMM::saveModel(string fileName){
	ofstream file(fileName);
	//file << "# N : Number of states" << endl;
	//file << N << endl;
	//file << "# M : Number of possible observations OR size of codebook" << endl;
	//file << M << endl;
	//file << "# The initial state distribution ( N x 1 Matrix)" << endl;
	printPI(file);
	//file << "# The state transition probability matrix ( N x N )" << endl;
	printA(file);
    //file << "# Observation probability distribution ( N x M Matrix)" << endl;
    printB(file);
    file.close();
}

void HMM::setObservationSeq(vector <unsigned int> arr, int T){
	O.clear();
	for (int i = 0; i <= T; ++i)
    {
        O.push_back(0);
    }
    for (int i = 0; i < T; ++i)
    {
        O[i] = arr[i];
    }
}
/*
long double HMM::predict(vector <unsigned int> obs, int T)
{
	init(T);
	int i,j,t;
    for(i=0;i<N;i++){ 
        alpha[0][i] = pi[i] * B[i][obs[0]];
    }

    for(t=0; t<T-1; t++)
    {
        for(j=0; j<N; j++)
        {
            alpha[t+1][j] = 0;
            for(i=0; i<N; i++)
            {
                    alpha[t+1][j] += alpha[t][i] * A[i][j];
            }
			//cout<<obs[1];
			//cout<<obs[26<<endl;
            alpha[t+1][j] *= B[j][obs[t+1]];
        }
    }

    prob_observ_given_model = 0.0;
    for(i=0; i<N; i++){
        prob_observ_given_model += alpha[T-1][i];
    }
	
	return prob_observ_given_model;
}
*/
void HMM::calculateAlpha(int T){
    int i,j,t;
    for(i=0;i<N;i++){ 
        alpha[0][i] = pi[i] * B[i][O[0]];
    }

    for(t=0; t<T-1; t++)
    {
        for(j=0; j<N; j++)
        {
            alpha[t+1][j] = 0;
            for(i=0; i<N; i++)
            {
                    alpha[t+1][j] += alpha[t][i] * A[i][j];
            }
            alpha[t+1][j] *= B[j][O[t+1]];
        }
    }
	
    prob_observ_given_model = 0.0;
    for(i=0; i<N; i++){
        prob_observ_given_model += alpha[T-1][i];
    }
}

void HMM::calculateBeta(int T){
    for (int i = 0; i < N; ++i)
    {
        beta[T-1][i] = 1;
    }
    
    for(int t = T-2; t>=0 ; t--)
    {   
        for (int i = 0; i < N; ++i)
        {
            beta[t][i] = 0;
            for (int j = 0; j < N; ++j)
            {
                beta[t][i] += beta[t+1][j]*A[i][j]*B[j][O[t+1]];
            }
        }
    }
}

void HMM::calculateGamma(int T){

    for(int t=0; t<T ;t++)
    {
        type_prob denom = 0;
        for(int i=0; i<N ;i++)
        {
			gamma[t][i]=0;
			for(int j=0;j<N;j++)
			{
				gamma[t][i] += xi[t][i][j];
			}            //denom += alpha[t][i] * beta[t][i];
        }
        
    }
}

void HMM::viterbi(int T){
    // Initialise
	
    for (int i = 0; i <N ; ++i)
    {
        delta[0][i] = pi[i]*B[i][O[0]];
        psi[0][i] = 0;
    }

    for (int j = 0; j <N ; ++j)
    {
        for (int t = 1; t <T ; ++t)
        {
            delta[t][j] = 0;
            for (int i = 0; i <N-1 ; ++i)
            {
                if( (delta[t][j] - delta[t-1][i]*A[i][j]) < ERROR){
                    delta[t][j] = delta[t-1][i]*A[i][j];
                    psi[t][j] = i;
                }
            }
            delta[t][j] *= B[j][O[t]];
        }
    }
    pStar = 0;
    for (int i = 0; i < N; ++i)
    {
        if( (pStar - delta[T-1][i]) < ERROR){
            pStar = delta[T-1][i];
            qStar[T-1] = i;
        }
    }

    for (int t = T-2; t < 1; t--)
    {
        qStar[t] = psi[t+1][ qStar[t+1] ];
    }
}

void HMM::calculateXi(int T){
    for(int t = 0; t <T-1 ; t++){
        type_prob denom = 0;
        type_prob expression;
        for (int i = 0 ;i < N ; ++i)
        {
            for (int j = 0; j < N; ++j)
            {
                expression = alpha[t][i]*A[i][j]*B[j][O[t+1]]*beta[t+1][j];
                xi[t][i][j] = expression;
                denom += expression;
            }
        }

        for (int i = 0; i < N ; ++i)
        {
            for (int j = 0; j < N; ++j)
            {
               if (denom!=0)
				   xi[t][i][j] /= denom;
			   else
				   xi[t][i][j]=0;
            }
        }
    }
}

void HMM::updateModel(int T){
 for (int i = 0; i < N; i ++)
    {
        pi[i] = gamma[0][i];

        for (int j = 0; j < N; j ++)
        {
            type_prob num_A = 0;
            type_prob denom_A = 0;
            for (int t = 0; t < T - 1; t ++)
            {
                    num_A += xi[t][i][j];
                    denom_A += gamma[t][i];
            }
            A[i][j] = num_A / denom_A;
            //if(A[i][j] == 0 )
               // A[i][j] = MIN_VALUE;
        }
        
        for(int k=0; k< M; k++)
        {
            type_prob num_B = 0;
            type_prob denom_B = 0;
            for(int t=0; t<T; t++)
            {
                    if(O[t] == k ) num_B += gamma[t][i];
                    denom_B += gamma[t][i];
            }
			
            B[i][k] = num_B / denom_B;
        }
    }
}

template<class In>
void HMM::printAlpha(In& file, int T){
    for (int t = 0; t <T ; ++t)
    {
        for (int i = 0; i <N ; ++i)
        {
            file << alpha[t][i] << '\t';
        }
        file << endl;
    }
    file << endl;
}

template<class In>
void HMM::printBeta(In& file, int T){
    for (int t = 0; t <T ; ++t)
    {
        for (int i = 0; i <N ; ++i)
        {
            file << beta[t][i] << '\t';
        }
        file << endl;
    }
    file << endl;
}

template<class In>
void HMM::printGamma(In& file, int T){
    for (int t = 0; t <T ; ++t)
    {
        for (int i = 0; i <N ; ++i)
        {
            file << gamma[t][i] << '\t';
        }
        file << endl;
    }
    file << endl;
}

template<class In>
void HMM::printXi(In& file, int T){
	//xi was t-1
    for (int t = 0; t <T ; ++t)
    {
        file << "t = " << t << endl;
        for (int i = 0; i <N ; ++i)
        {
            for (int j = 0; j <N ; ++j){
                file << xi[t][i][j] << '\t';
            }
        }
        file << endl;
    }
    file << endl;
}

template<class In>
void HMM::printPsi(In& file, int T){
    for (int t = 0; t <T ; ++t)
    {
        for (int i = 0; i <N ; ++i)
        {
            file << psi[t][i] << '\t';
        }
        file << endl;
    }
    file << endl;
}

template<class In>
void HMM::printDelta(In& file, int T){
    for (int t = 0; t <T ; ++t)
    {
        for (int i = 0; i <N ; ++i)
        {
            file << delta[t][i] << '\t';
        }
        file << endl;
    }
    file << endl;
}

template<class In>
void HMM::printA(In& file){
    for(int i=0 ; i<N ; i++)
    {
        for (int j = 0; j <N ; ++j)
        {
            file << A[i][j] << '\t';
        }
        file << endl;
    }
    file << endl;
}

template<class In>
void HMM::printB(In& file){
    for(int i=0 ; i<N ; i++)
    {
        for (int j = 0; j <M ; ++j)
        {
            file << B[i][j] << '\t';
        }
        file << endl;
    }
    file << endl;
}

template<class In>
void HMM::printPI(In& file){
    for(int i=0 ; i<N ; i++){
        file << pi[i] << '\t';
    }
    file << endl;
}

template<class In>
void HMM::printALL(In& file, int T){
	file << "# A " << endl;
    printA(file);
    file << "# B " << endl;
    printB(file);
    file << "# Pi " << endl;
    printPI(file);
    file << "# Alpha " << endl;
    printAlpha(file,T);
    file << "# Beta " << endl;
    printBeta(file,T);
    file << "# Gamma " << endl;
    printGamma(file,T);
    file << "# Xi " << endl;
    printXi(file,T);
    file << "# Psi " << endl;
    printPsi(file,T);
    file << "# Delta " << endl;
    printDelta(file,T);
}

void HMM::train(vector <unsigned int> arr, int T){
	init(T);
    setObservationSeq(arr,T);
    int i = 0;
    type_prob prevProb = 0;
    do{
		if(i>1000)
			break;
		ofstream file("Logs\\iter_"+ toString(i) +".txt" , ios::out);
        if (!file.is_open()){
            cout << "Cannot open File ";
            exit(0);
        }
		init(T);
        viterbi(T);
        prevProb = prob_observ_given_model;
		calculateAlpha(T);
        calculateBeta(T);
        calculateXi(T);
		calculateGamma(T);
        updateModel(T);
        cout <<"prob : "<< prob_observ_given_model << " Iteration : "<< i ;
        cout <<" Diff "<< prevProb - prob_observ_given_model << endl;
        i++;
        printALL(file, T);
        file.close();
    }while	(abs(prevProb -prob_observ_given_model ) > ERROR || i<10);
	//while(i<20);
		//
	// while(i<20);
}

void train()
{
	//ShellExecute(0, _T("open"), _T("exec.exe"), _T("5 record.wav record.txt"), 0, SW_SHOWNORMAL);
	//Sleep(5000);
	//cout << "executing" << endl;
	//string filename = "record.txt";
	
	//cout<<"The observation sequence is being taken from data\\assgn\\*.txt"<<endl;
	

	double max = 0;
	int word=0, t=0;
	for(int fin=0;fin<10;fin++)
	{
		vector< vector<type_prob>> Aavg;
		Aavg.resize(N+1);
		for(int i=0;i<N;i++)
		{
			Aavg[i].resize(M+1);
		}
		vector< vector<type_prob>> Bavg;
		Bavg.resize(N+1);
		for(int i=0;i<N;i++)
		{
			Bavg[i].resize(M+1);
		}
		vector <type_prob> piavg;
		piavg.resize(N+1);
		
		for(int p=0;p<N;p++)
		{
			piavg[p]=0;
			for(int q=0;q<M;q++)
			{
				Aavg[p][q]=0;
				Bavg[p][q]=0;
			}
		}
		for (int i=0;i<10;i++)
		{	
			cout<<"For "<<fin<<": Model "<<i<<":"<<endl;
			ifstream readob("data\\quantized\\"+toString(fin)+"\\"+toString(i+1)+".txt",ios::in);
			unsigned int tem=0;
			vector <unsigned int> observationseq;
			int count=0;
			while(readob>>tem)
			{
				observationseq.push_back(tem);
				//cout<<tem;
				count++;
			}
			componentmodels[i]=new HMM();
			//componentmodels[i]->saveModel("models\\"+toString(fin)+" "+toString(i+1)+"_initialcompmodel.txt");
			componentmodels[i]->train(observationseq, observationseq.size());	
			//componentmodels[i]->saveModel("outputs\\"+toString(i+1)+"finalmodel.txt");
			cout<<endl;
			for(int p=0;p<N;p++)
			{
				piavg[p]+=pi[p];
				for(int q=0;q<M;q++)
				{
					Aavg[p][q]+=A[p][q];
					Bavg[p][q]+=B[p][q];
				}
			}
		}

		
		
		for(int p=0;p<N;p++)
		{
			pi[p]=piavg[p]/10;
			for(int q=0;q<M;q++)
			{
				A[p][q]=Aavg[p][q]/10;
				B[p][q]=Bavg[p][q]/10;
			}
		}
		
		finalmodels[fin]->saveModel("models\\"+toString(fin)+".txt");
		//finalmodels[fin] = new HMM(Aavg, Bavg, piavg);
	}

	
	
}

/*void loadModel(string filepath)
{
	ifstream model(filepath.c_str(),ios::in);
	HMM* mo = new HMM();
	//mo->init(
	string line;
	int which=0;
	if(model.is_open())
	{
		cout<<"open"<<endl;

	while(getline(model,line))
	{
		cout<<line<<endl;
		if(line[0]=='#')
		{
			which++;
			if(which==1 || which==2)
			{
				continue;
			}
			else if(which==3)
			{
				for(int i=0;i<N;i++)
				{
					model>>pi[i];
				}
			}
			else if(which==4)
			{
				
				for (int row = 0; row < N; row++)
				{
					for (int col = 0; col < M; col++)
					{   
						model >> A[row][col];
					 }   
				}
			}
			else if (which==5)
			{
				for (int row = 0; row < N; row++)
				{
					for (int col = 0; col < M; col++)
					{   
						model >> B[row][col];
					 }   
				}
			}

			continue;
		}
		else if(which==1)
		{
			N=stoi(line);
		}
		else if(which==2)
		{
			M=stoi(line);
		}
	}

	}
}*/

void LoadModel(int k, int j, int T)
{
	ifstream inp;
	string filename = "models\\"+toString(k)+".txt";
	inp.open(filename.c_str());
	A.clear();
	B.clear();
	O.clear();
	pi.clear();
	finalmodels[(k+1)*j]=new HMM();
	finalmodels[(k+1)*j]->init(T);
	
		long double temp;
		for(int i=0; i<N; i++)
		{	
			inp >> temp;
			pi[i]=temp;
		}
		for(int i=0; i<N; i++)
		{
			for(int j=0; j<N; j++)
			{
				inp >> temp;
				A[i][j]=temp;
			}
		}
		for(int i=0; i<N; i++)
		{
			for(int j=0; j<M; j++)
			{
				inp >> temp;
				B[i][j]=temp;
			}
		}
	
	return;
}
void test(){
	cout<<"Models saved in models folder"<<endl;
	
	int correct = 0;
		for(int i=0;i<=9;i++)
			{
				for(int j=1;j<=5;j++)
				{	ifstream readob("test\\"+toString(i)+toString(j)+".txt",ios::in);
					unsigned int tem;
					vector <unsigned int> observationseq;

					while(readob>>tem)
					{
						observationseq.push_back(tem);
						//cout<<tem;
					}
					long double prob=0;
					long double maxprob = INT_MIN;
					readob.close();
					//int maxmodel;
					cout<<endl;
					int maxword = -1;
					for(int k=0;k<10;k++)
					{
						LoadModel(k,j,observationseq.size());

						finalmodels[((k+1)*j)]->init(observationseq.size());
						finalmodels[((k+1)*j)]->setObservationSeq(observationseq,observationseq.size());
						finalmodels[((k+1)*j)]->calculateAlpha(observationseq.size());
						prob = prob_observ_given_model;
						cout<<"Prob of word being "<<k<<"is: "<<prob<<endl; 
						if(prob>maxprob)
						{
							maxprob= prob;
							maxword = k;
						}
					}
					cout<<"Prediction for the word is "<<maxword<<" with prob: "<<maxprob<<endl;
					if(i==maxword)
					{
						correct++;
					}
				}
		}
	
	
}

int main()
{
	train(); 
	//test();
	//loadModel("outputs\\14_initialcompmodel.txt");
	cout<<"done"<<endl;
	int a; cin>>a;
	return 0;
}