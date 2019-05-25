#ifndef GPAgraph_H
#define GPAgraph_H

#include "Graph.h"
#include "Helperfunctions.h"
#include <vector>
#include <algorithm>
#include <cmath>
#include <random>
#include <fstream>


#define PI 3.1415926

double getMinBucketDistEstimate(std::vector<std::vector<std::vector<std::vector<double> > > > &bucketSarkok,
	std::vector<double> p, int i, int j, double unit);


double getMaxBucketDist(std::vector<std::vector<std::vector<std::vector<double> > > > &bucketSarkok, std::vector<double> p, int i, int j);


// It calculates the expected number of nodes with degree k after inserting n nodes to the PA model 
double getProb(int n, int m, double delta, int k);



// calculates the value of F function, in the case when F = F0
double getFvalue0(double r, double u);


// calculates the value of F function, in the case when F = F1
double getFvalue1(double r, double u);



// calculates the value of F function, in the case when F = F2
double getFvalue2(double beta, double u, int n, double theta);


// calculates the value of F function, in the case when F = F3
double getFvalue3(double beta, double u);



// calculate the integral which appered in the GPA model, and was denoted by I
double getI(double beta, int n);



// generating graph from the model when F = F1
// in this case we use S = [0,1]^2 and the torus metric
// it is really similar to the case when we are in the sphere
// There is a very effective way to generate such graphs

// A brute force solution can be that in each step you generate a point on the sphere,
// after you calculate the distance from all existing nodes, then store the ones which are 
// close enough, and after generate the edges in the correct way

// IMPROVEMENT: Instead of calculating the distance from each point in each step, we can make 
// buckets, and we can store from each nodes which bucket they belong.
// We can calculate the distances between the buckets, intead of each nodes.
// There are 3 possibilities for each bucket:
//  1. if the whole bucket is outside of the given radius, then we can ignore those buckets
//  2. if the whole bucket is inside of the given radius, then we know that the nodes in the bucket are in the radius
//  3. if the bucket has parts in and out as well, then we have to check one by one each node in the bucket whether they are close enough.

// Generating this way allows us to get graphs with more then 10^6 nodes very fast

Graph generateGPA1(int n, int m, double delta, double r, int k, std::mt19937& gen);

// to generate PA model we use the fact, that in the GPA model with F=F1, if r>sqrt(2)/2 then we get the simple PA model
Graph generatePA(int n, int m, double delta,  std::mt19937& gen);


// generating graph from the model when F = F2
// this is on the sphere, with the angular metric
Graph generateGPA2(int n, int m, double alpha, double delta, double beta, std::mt19937& gen);


// generating graph from the model when F = F3
// this is on the sphere, with the angular metric
Graph generateGPA3(int n, int m, double alpha, double delta, double beta, std::mt19937& gen);



// simulate PA graphs with the given parameters
void simulatePAandSave(std::string path, int  T, std::vector<int> N, std::vector<int> M, std::vector<double> deltas, std::mt19937& gen);



// simulate GPA graphs with the given parameters, F=F1
void simulateGPA1andSave(std::string path, int  T, std::vector<int> N, std::vector<int> M, std::vector<double> R,
 std::vector<double> deltas, std::mt19937& gen);


// simulate GPA graphs with the given parameters, F=F2
void simulateGPA2andSave(std::string path, int  T, std::vector<int> N, std::vector<int> M, std::vector<double> betas,
	std::vector<double> alphas, std::vector<double> deltas, std::mt19937& gen);


// simulate GPA graphs with the given parameters, F = F3
void simulateGPA3andSave(std::string path, int  T, std::vector<int> N, std::vector<int> M, std::vector<double> betas,
	std::vector<double> alphas, std::vector<double> deltas, std::mt19937& gen);







#endif
