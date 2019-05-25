#include <iostream>
#include <vector>
#include <map>
#include <time.h>
#include <fstream>

#include "Graph.h"
#include "GPAgraph.h"
#include "Helperfunctions.h"


using namespace std;

#define PI 3.1415926

int main()
{
	srand((unsigned)time(0));
	clock_t start = clock();
	std::mt19937 gen;
	gen.seed(time(0));//if you want different results from different runs

	string path = "E:/egyetem/msc szakdoga/szimulaciok/GPA2_proba";


	int T = 10; // how many simulation do you want for each 
	//parameters
	vector<int> N{ 50 }; 
	vector<int> M{ 5 };
	vector<double> betas{ 0.5,1,1.5,2,3 };
	vector<double> alphas{ 0 };
	vector<double> deltas{ 0 };

	simulateGPA2andSave(path, T, N, M, betas, alphas, deltas,gen);


	
}