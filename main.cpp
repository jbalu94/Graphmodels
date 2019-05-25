#include <iostream>
#include <vector>
#include <map>
#include <time.h>
#include <fstream>

#include "Graph.h"
#include "GPAgraph.h"
#include "HPAgraph.h"
#include "Helperfunctions.h"


using namespace std;

#define PI 3.1415926

int main()
{
	srand((unsigned)time(0));
	clock_t start = clock();
	std::mt19937 gen;
	gen.seed(time(0));//if you want different results from different runs

	// add the path to your folder where you want to store the results
	string path = "E:/egyetem/msc szakdoga/szimulaciok";


	int T = 10; // how many simulation do you want for each 

	//set the parameters
	vector<int> N{ 500 }; 
	vector<int> M{ 5 };
	vector<double> betas{ 0.5 };
	vector<double> alphas{ 0 };
	vector<double> deltas{ 0.1,1,10 };
	vector<double> R{ 0.1,0.2,0.3,0.4,0.5,0.6 };


	// choose the model you wish to run

	//simulatePAandSave(path, T, N, M, deltas, gen);
	//simulateGPA1andSave(path, T, N, M, R, deltas, gen);
	//simulateGPA2andSave(path, T, N, M, betas, alphas, deltas,gen);
	//simulateGPA3andSave(path, T, N, M, betas, alphas, deltas, gen);
	//simulateSPAandSave(path, T, N, M, betas, alphas, deltas, gen);
	simulateHPAandSave(path, T, N, M, betas, deltas, gen);
	//Graph G = generateHPA(20, 3, 0.66, 1, gen);

	return 0;

}