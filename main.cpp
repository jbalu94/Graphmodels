#include <iostream>
#include <vector>
#include <fstream>
#include <time.h>

#include "Graph.h"
#include "GPAgraph.h"
#include "HPAgraph.h"
#include "Helperfunctions.h"


int main()
{
	srand((unsigned)time(0));
	std::mt19937 gen;
	gen.seed(time(0));//if you want different results from different runs

	// add the path to your folder where you want to store the results
	std::string path = "E:/egyetem/msc szakdoga/szimulaciok";


	int T = 10; // how many simulation do you want for each 

	//set the parameters
	std::vector<int> N{ 2000 }; 
	std::vector<int> M{ 5 };
	std::vector<double> betas{ 0.2,0.4,0.6,0.8 };
	std::vector<double> alphas{ 0 };
	std::vector<double> deltas{-4,-2,0,4,10};
	std::vector<double> R{ 0.1,0.2,0.3,0.4,0.5,0.6 };


	//------------
	// choose the model you wish to run
	//------------

	simulatePAandSave(path, T, N, M, deltas, gen);
	//simulateGPA1andSave(path, T, N, M, R, deltas, gen);
	//simulateGPA2andSave(path, T, N, M, betas, alphas, deltas,gen);
	//simulateGPA3andSave(path, T, N, M, betas, alphas, deltas, gen);
	//simulateSPAandSave(path, T, N, M, betas, alphas, deltas, gen);
	//simulateHPAandSave(path, T, N, M, betas, deltas, gen);



	//-------------
	// in case you would like to generate one graph with certain paramters, and save it
	//-------------
	//Graph G = generateHPA(500, 3, 0.75, 1, gen);
	//writeToFile(G, path+"/HPAgraph_500_3_0.75_1.txt");

	return 0;

}