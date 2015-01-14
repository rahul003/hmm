// Assignment3.cpp : Defines the entry point for the console application.
#include <iostream>
#include <stdio.h>
#include <string>
#include <fstream>
#include <time.h>
#include <windows.h>
#include <vector>
#include <sstream>
//Globals
using namespace std;
//ifstream inp;
ofstream cep;//("data\\0\\cep.txt",fstream::app);
	
int totsam;
string fname;
double *variance;
int point=0;
int highest =0;
vector <string> filelist;
int windowSize;
double normMax;
int k;
double minError;
	
string toStringh(int number)
{
   stringstream ss;//create a stringstream
   ss << number;//add number to the stream
   return ss.str();//return a string with the contents of the stream
}
void GetAllFiles(string dirname)
{
	HANDLE hFind;
	WIN32_FIND_DATAA data;
	dirname = "data\\"+dirname;
	hFind = FindFirstFileA(dirname.c_str(), &data);
	if (hFind != INVALID_HANDLE_VALUE) {
	  do {
		//printf("%s\n", data.cFileName);
		filelist.push_back(data.cFileName);
	  } while (FindNextFileA(hFind, &data));
	  FindClose(hFind);
	}

}
double* normalize(double* arr, int count, double max, double min, double normMax)
{
	double modmax;
	if (max > -min)
		modmax = max;
	else
		modmax = -min;
	double ratio = normMax / modmax;
	for (int i = 0; i<count; i++)
	{
		arr[i] = arr[i] * ratio;
	}
	//std::cout<<"Normalization carried out"<<endl;
	return arr;
}
double* shift(double* arr, double sum, int count)
{
	double shift = sum / count;
	//std::cout << "Average of all samples is " << shift << endl;
	if ((shift>0 && shift < 1) ||( shift<0 && shift >-1))
	{
		//std::cout << "Dc shift is not much, so not shifting waveform" << endl;
		
	}
	//shifting
	else
	{
		//double sum = 0;
		for (int i = 0; i<count; i++)
		{
			arr[i] -= shift;
			//sum = sum + arr[i];
		}
		//cout << sum / count << endl;
	}


	return arr;
}
void getRs(double* arr, int start, int end, int p, double* rval)
{	//get the autocorrelation Ris 
	//cout << "ris" << endl;
	for (int k=0;k<=p;k++)
	{
			rval[k]=0;
			for(int m=start; m<end-k; m++)
			{

				rval[k] +=arr[m]*arr[m+k];
				//if (start == 3600)
					//cout << rval[k] << endl;
			}
			//rval[k]/=(end-k-start+1);
			//cout << rval[k] << endl;
	}
	
	return ;
}
double sumoverm(double* a, double* c, int m)
{	//implement the sum over m term while calculating Cis
	int i=0;
	double sum=0;
	for(i=1;i<m;i++)
	{sum=sum + (i/(double)m)*c[i]*a[m-i];
	}
	//cout << "sum " << sum << endl;
	return sum;
}
void getCis(double* a, int p, double* c)
{	//implement the equation which converts Ai to Ci
	int i=0;
	for(i=1;i<=12;i++)
	{
		c[i]=a[i] + sumoverm(a, c, i);
		//cout << c[i] << endl;
	}

	return ;
}
void getAis(double* R, int p, double * a)
{
	double alpha[13][13];
	for(int i =0;i<=12;i++)
	{
		for(int j=0;j<=12;j++)
		{
			alpha[i][j]=0;
		}
	}

	double E[13];
	double k[13];
	E[0] = R[0];
	for(int i=1; i<=p; i++)
	{
		double sum =0;
		for(int j=1; j<i; j++)
		{
			sum += alpha[j][i-1]*R[i-j];
		}
		
		k[i] = (R[i] - sum)/E[i-1];
		E[i] = (1 - k[i]*k[i])*E[i-1];
		alpha[i][i] = k[i];
		for(int j=1; j<i; j++)
		{
			alpha[j][i] = alpha[j][i-1] - k[i]*(alpha[i-j][i-1]);
		}
		
	}
	for(int i=1; i<=p; i++)
	{
		a[i] = alpha[i][p];
		//cout << a[i - 1] << endl;
	}
	return ;
}
double* getCepstrals(double* arr, int start, int end)
{
	//handles all subroutines and prints Cis ultimately
	double* R = new double[13];
	getRs(arr,start, end, 12, R);
	double* a = new double[13];
	getAis(R, 12,a);

	double* c = new double[13];
	getCis(a, 12, c);
	/*
	fstream outp ("outputs\\_out.txt", fstream::app);
	for(int i=1; i<=12; i++)
	{
		outp << *(c+i) << "," ;
	}
	outp << endl; */
	return c;
}
double calculate_tokura(double *vector1, double *vector2)
{
double dist = 0;
	for (int i = 0; i <12; i++)
		dist += ((vector1[i] - vector2[i])*(vector1[i] - vector2[i]));/// variance[i];// 
	return dist;
}
void findVariance(double** cepstral_matrix, double* variance_val, int num_frames)
{
	double average[12]={0};
	int i, j;
	for (i = 0; i < 12; i++)
	{
		for (j = 0; j < num_frames; j++)
		{
			average[i] += cepstral_matrix[j][i];
		}
		average[i] /= num_frames;
	}

	for (i = 0; i < 12; i++)
	{
		for (j = 0; j < num_frames; j++)
		{
			variance_val[i] += (cepstral_matrix[j][i] - average[i])*(cepstral_matrix[j][i] - average[i]);
		}
		variance_val[i] /= num_frames;
		//cout<<"var:"<<variance_val[i]<<endl;
	}

	return;

}
double** convert2D(double* cepstral_data, int num_frames)
{
	double **data2D = new double*[num_frames];
	for (int i = 0; i < num_frames; ++i)
		data2D[i] = new double[12];
	for (int i = 0; i < num_frames; i++)
	{
		for (int j = 0; j < 12; j++)
		{
			data2D[i][j] = cepstral_data[(i * 12) + j];
		}
	}

	return data2D;

}
void assign_cluster(double **data, double** clusters, int num_frames, int k, int** cluster_ele)
{
	int index=0; 
	//int **cluster_ele = new int*[k];
	//for (int i = 0; i < k; i++)
		//cluster_ele[i] = new int[num_frames];
	for (int i=0;i<k;i++)
	{
		for(int j=0;j<num_frames;j++)
		{
			cluster_ele[i][j]=0;
		}
	}
	for (int i = 0; i < num_frames; i++)
	{
		double min = INT_MAX, temp;
		index = 0;
		for(int j = 0;  j< k; j++)
		{
			temp = calculate_tokura(data[i], clusters[j]);
			if (min > temp)
			{
				min= temp;
				index= j;
			}
		}
		cluster_ele[index][i] = 1;	
	}
	return;
}
double calc_error(double ** centroids, double ** cepstral_data, int ** cluster_indices, int k, int num_frames, double * variance)
{
	
	double error = 0, clusterError=0, maxErr=0, temp;
	int val=0;
	//cout<<"k is"<<k<<num_frames<<endl;
	for (int cluster = 0; cluster < k; cluster++)
	{
		clusterError = 0;
		for(int j=0;j<num_frames;j++)
		{
			if(cluster_indices[cluster][j])
			{
				temp=calculate_tokura(centroids[cluster], cepstral_data[j]);
				error+= temp;
				clusterError+=temp;
				val= j;
			}
		}
		if (clusterError>maxErr)
		{
			maxErr=clusterError; 
			point= val; 
			highest=cluster;
		}
	}
	return error;
}
void new_centroid(double ** centroids, double ** cepstral_data, int** cluster_indices, int num_frames, int k)
{
	double* sum = new double[k * 12];
	int* clustersize = new int [k];
	for(int temp=0;temp<k;temp++)
	{
		clustersize[temp]=0;
	}
	//memset(clustersize, 0, k*sizeof(double));
	for (int cluster = 0; cluster < k; cluster++)
	{
		memset(sum, 0, 12*k*sizeof(double));
		for(int j=0;j<num_frames;j++)
		{
			if(cluster_indices[cluster][j])
			{
				for(int dim=0;dim<12;dim++)
				{
					sum[cluster * 12 + dim] += cepstral_data[j][dim];
				}
				clustersize[cluster]++;
					
			}
		}
		
/*
		while (cluster_indices[cluster][clustersize[cluster]] != -1)
		{
			for (int dim = 0; dim < 12; dim++)
			{
				sum[cluster * 12 + dim] += cepstral_data[cluster_indices[cluster][k]][dim];
			}
			clustersize[cluster]++;
		}*/
		for (int dim = 0; dim < 12; dim++)
		{
			if(clustersize[cluster]==0)
			{
				//cout<<"reassign"<<endl;
			   centroids[cluster][dim]= cepstral_data[point][dim];
			   cluster_indices[cluster][point]=1;
			   cluster_indices[highest][point]=0;
			}
			
			else
				{
					centroids[cluster][dim] = sum[cluster * 12 + dim] / clustersize[cluster];
					//cout<<centroids[cluster][dim]<<" ";
			}
		}
		//cout<<endl;
	}
	return;
}
void K_Means(int k, double* data, int num_frames, double min_error, int digit,int fileid)
{
	double **cepstral_data = convert2D(data, num_frames);

	//choose k-random points
	//variance = new double[12];
	//findVariance(cepstral_data, variance, num_frames);
	int i, j;
	double **centroids = new double*[k];
	for (int i = 0; i < k; ++i)
		centroids[i] = new double[12];
	srand( (unsigned)time( NULL ) );
	for (i = 0; i < k; i++)
	{
		j = rand() % num_frames;
		for (int dim = 0; dim < 12; dim++)
		{
			centroids[i][dim] = cepstral_data[j][dim];
		}
	}
	double error=INT_MAX;

	int **cluster_indices = new int*[k];
	for (int i = 0; i < k; i++)
		cluster_indices[i] = new int[num_frames];
	
	fstream outp("outputs\\kmeans"+toStringh(digit)+","+toStringh(fileid)+"_errors.txt", fstream::app);
	fstream log("logs\\kmeans"+toStringh(digit)+","+toStringh(fileid)+"_log.txt", fstream::app);
	assign_cluster(cepstral_data, centroids, num_frames ,k, cluster_indices);
	error = calc_error(centroids, cepstral_data, cluster_indices, k, num_frames, variance);
	double preverror = error;
	int num_iter=0;
	do
	{
		new_centroid(centroids, cepstral_data, cluster_indices,num_frames, k);
		assign_cluster(cepstral_data, centroids, num_frames ,k,cluster_indices);
		for(int i =0;i<k;i++)
		{
			int num=0;
			for(int j=0;j<num_frames;j++)
			{
				if(cluster_indices[i][j])
					num++;
			}
			log<<"num of points to cluster "<<i<<" is "<<num<<endl;
		}
		preverror =error;
		error = calc_error(centroids, cepstral_data, cluster_indices, k, num_frames, variance);
		outp << error/num_frames << endl;
		log<<"error "<<error/num_frames<<endl;
		//cout<<"PrevError: "<<preverror<<endl;
		cout<<num_iter<<", Error: "<<(error/num_frames)<<endl;
		num_iter++;
	//}while(num_iter<100);
	}while (abs(preverror-error)>min_error);

	//hardcoded num of class
	int total=0;
	int current=0;

	ofstream assignment("data\\quantized\\"+toStringh(digit)+"\\"+toStringh(fileid)+".txt", fstream::app);
	for(int i=0;i<num_frames;i++)
	{
		for(int j=0;j<k;j++)
		{
			if(cluster_indices[j][i])
				assignment<<j<<endl;
		}

	}
	assignment.close();
	log.close();
}
void makeCepstral(string path)
	{
	
	//GetAllFiles(dir);

	//		string num = "";
		//num=num+i;
		//string path =parpath+"\\"+toStringh(i)+".txt"; 
		ifstream inp(path.c_str(), fstream::in);
		double sample;
		double sum = 0;
		int count = 0;
		if(inp.is_open())
		{
			//cout<<"entered"<<endl;
			while (inp >> sample)
			{
				sum = sum + sample;
				count++; //count no. of samples above
			}

			inp.clear();
			inp.seekg(0, ios::beg); //reset filepointer to beginning of file
		
			double* arr = new double[count]; //dynamically allocate array
			int i = 0;
			double max = 0;
			double min = 0;
			while (inp >> sample)
				{
					//store values to array
					arr[i] = sample;
					if ((sample > 0) && (sample > max))
						max = sample;
					if ((sample < 0) && (sample < min))
						min = sample;
					i++;
				}
			arr = shift(arr, sum, count);
			arr = normalize(arr, count, max, min, normMax);
		
			//cout<<"count "<<count<<endl;
			double * train_set = new double [(count / windowSize) * 4 *12];
			int num = 0;
			double * temp;
			int co=0;

			for (int i = 0; i < count-windowSize; i+=(windowSize/4))
			{
				co++;
				//cout << "window " << i << " " << i + windowSize<<endl;
				temp = getCepstrals(arr, i, i + windowSize);
				for (int j = 1; j <= 12; j++)
				{
					//train_set[num*12+j] = temp[j];
					//cout<<temp[j]<<endl;
					cep << temp[j] << endl;
				}
				num++;
			}
		//cout<< "numframes"<<co<<endl;
		}

		inp.close();
	
	
}
/*
LBG Code follows
*/

double ** new_centroid_LBG(double ** centroids, double ** cepstral_data, int** cluster_indices, int num_frames, int k)
{
	double* sum = new double[k * 12];
	int* clustersize = new int [k];
	for(int temp=0;temp<k;temp++)
	{
		clustersize[temp]=0;
	}
	for (int cluster = 0; cluster < k; cluster++)
	{
		memset(sum, 0, 12 * k*sizeof(double));
		for (int j = 0; j<num_frames; j++)
		{
			if (cluster_indices[cluster][j])
			{
				for (int dim = 0; dim<12; dim++)
				{
					sum[cluster * 12 + dim] += cepstral_data[j][dim];
				}
				clustersize[cluster]++;

			}
		}

		/*
		while (cluster_indices[cluster][clustersize[cluster]] != -1)
		{
		for (int dim = 0; dim < 12; dim++)
		{
		sum[cluster * 12 + dim] += cepstral_data[cluster_indices[cluster][k]][dim];
		}
		clustersize[cluster]++;
		}*/
		for (int dim = 0; dim < 12; dim++)
		{
			if (clustersize[cluster] == 0)
			{
				//cout<<"reassign"<<endl;
				centroids[cluster][dim] = cepstral_data[point][dim];
				cluster_indices[cluster][point] = 1;
				cluster_indices[highest][point] = 0;
			}

			else
			{
				centroids[cluster][dim] = sum[cluster * 12 + dim] / clustersize[cluster];
				//cout<<centroids[cluster][dim]<<" ";
			}
		}
		//cout<<endl;
	}
	return centroids;
}
void get_centroid(double* cluster, double** cepstral_data, int num_frames)
{ 
	for (int dim = 0; dim<12; dim++)
		{
			cluster[dim] =0;
		}
	for (int j = 0; j<num_frames; j++)
	{
		for (int dim = 0; dim<12; dim++)
		{
			cluster[dim] += cepstral_data[j][dim];
		}
	}

	for (int dim = 0; dim < 12; dim++)
	{
		cluster[dim] = cluster[dim] / num_frames;
	}
}
int** assign_cluster_LBG(double **data, double** clusters, int num_frames, int k)
{
	int index = 0;
	int **cluster_ele = new int*[k];
	for (int i = 0; i < k; i++)
		cluster_ele[i] = new int[num_frames];
	for (int i = 0; i<k; i++)
	{
		for (int j = 0; j<num_frames; j++)
		{
			cluster_ele[i][j] = 0;
		}
	}
	for (int i = 0; i < num_frames; i++)
	{
		double min = INT_MAX, temp;
		index = 0;
		for (int j = 0; j< k; j++)
		{
			temp = calculate_tokura(data[i], clusters[j]);
			if (min > temp)
			{
				min = temp;
				index = j;
			}
		}
		cluster_ele[index][i] = 1;
	}
	return cluster_ele;
}
void LBG(int k, double* data, int num_frames, double min_error, int epsilon)
{
	double **cepstral_data = convert2D(data, num_frames);
	//choose k-random points

	variance = new double[12];
	findVariance(cepstral_data, variance, num_frames);
	fstream outp("outputs\\LBG_errors.txt", fstream::app);

	int i, j, n=1;
	double cluster[12], val;
	int** cluster_indices;
	 
	get_centroid(cluster, cepstral_data, num_frames);
	double **centroids = new double*[k];
	for ( i = 0; i < k; ++i)
		centroids[i] = new double[12];
	for ( i = 0; i < 12; i++)centroids[0][i] = cluster[i];

	fstream log("outputs\\lbg_log.txt", fstream::app);
	cluster_indices = assign_cluster_LBG(cepstral_data, centroids, num_frames, n);
	double error = calc_error(centroids, cepstral_data, cluster_indices, n, num_frames, variance);
	outp <<error << endl;
	while (n < k)
	{

		int next = n;
		for (i = 0; i < n; i++)
		{
			next++;
			for (j = 0; j < 12; j++)
			{
				if(next<k)
					{
						val = centroids[i][j];
						centroids[i][j] = val * (1 - epsilon);
						centroids[next-1][j] = val * (1+epsilon);
				}
			}
			
		}
		//n - current number of clusters
		n = next;
		cluster_indices = assign_cluster_LBG(cepstral_data, centroids, num_frames, n);

		for(int i =0;i<n;i++)
		{
			int num=0;
			for(int j=0;j<num_frames;j++)
			{
				if(cluster_indices[i][j])
					num++;
			}
			log<<"num of points to cluster "<<i<<" is "<<num<<endl;
		}
		centroids = new_centroid_LBG(centroids, cepstral_data, cluster_indices, num_frames, n);
		//cluster_indices = assign_cluster_LBG(cepstral_data, centroids, num_frames, n);
		double error = calc_error(centroids, cepstral_data, cluster_indices, n, num_frames, variance);
		outp <<error << endl;
		cout<<"cluster number = "<<n<<", Error:"<<error<<endl;
	}

}
int mainp()
{
	//get params from config file
	ifstream fparams("kmeansconfig.txt");
	fparams>>windowSize;
	fparams>>normMax;
	fparams >> k;
	fparams >> minError;
	fparams.close();
   
	cout<<"Cepstral coefficients being stored"<<endl;
	
	for (int i=0;i<10;i++)
	{	
		for(int j=1;j<=10;j++)
		{
				cep.open("data\\cep\\"+toStringh(i)+"\\"+toStringh(j)+"_cep.txt", fstream::app);
				makeCepstral( "data\\"+toStringh(i)+"\\"+toStringh(j)+".txt");
				cep.close();
		}
	}	

	for(int i=0;i<10;i++)
	{
		for (int j=1;j<=10;j++)
		{	
			//cout<<"Reading cepstral coeff for "<<<<endl;
			string fileToOpen ="data\\cep\\"+toStringh(i)+"\\"+toStringh(j)+"_cep.txt"; 
			ifstream cepinput(fileToOpen);
			double coeff;
			int num_coeffs = 0;
			while(cepinput>>coeff)
			{
				num_coeffs++;
			}
			cepinput.close();
			double* cepstraldata = new double[num_coeffs];
			cepinput.open(fileToOpen);
			int temp=0;
			while(cepinput>>coeff)
			{
				cepstraldata[temp]=coeff;
				temp++;
			}
			cout<<"kmeans running from now"<<endl;
			K_Means(8,cepstraldata,num_coeffs/12,0.00001,i, j);
			cout<<endl;
	
		}
	}	
	//cout<<"lbg running from now"<<endl;
	//LBG(128, cepstraldata, num_coeffs/12, 0.0001, 0.0001);
		
	int a;
	cout<<"Done"<<endl;
	cin >> a;
	
	return 0;
}
