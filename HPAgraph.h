#ifndef HPAgraph_H
#define HPAgraph_H

#include "Graph.h"
#include "Helperfunctions.h"
#include <vector>
#include <algorithm>
#include <cmath>
#include <random>
#include <fstream>



// calculate the hiperbolic distance between two points
double getHiperDist(std::vector<double> a, std::vector<double> b, double beta);


// generate the ang of the new node
double generateHPAang(Graph &G, int t, double beta, double delta, std::mt19937& gen, int maxi);


// generate a HPA graph 
Graph generateHPA(int n, int m, double beta, double delta, std::mt19937& gen);


// simulate HPA graphs with certain parameters, and save it to files
void simulateHPAandSave(std::string path,int T,std::vector<int> N, std::vector<int> M, 
	std::vector<double> betas, std::vector<double> deltas, std::mt19937& gen);


#endif