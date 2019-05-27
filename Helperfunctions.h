#ifndef Helperfunctions_H
#define Helperfunctions_H

#include <vector>
#include "Graph.h"

// this file contains the function which calculates the attributions of the certain graphs


double min(double a, double b);
double max(double a, double b);

//distance in the torus metric 
double distTorus(std::vector<double> a, std::vector<double>  b);



// dot product of two vector
double dot(std::vector<double> a, std::vector<double> b);

// calculate the angular between two vector
double getAng(std::vector<double> a, std::vector<double> b);

// generate a random floating number between [0,1]
double fRand(double fMin, double fMax);


// rankify a vector
std::vector<double> rankify(std::vector<double> & X);


// generating a random point on a sphere
std::vector<double> generateRandomToSphere(double r);


// calculate degree distribution
std::vector<int> getDegreeDist(Graph &G);

// calculate global clustring coefficient
double getGlobalClustering(Graph &G);

// calculate avg clustering coefficient
double getAvgClustering(Graph &G);

// calculate C(d) function, which is the mean of the local clustering coefficient which has 
std::vector<double> getCd(Graph& G);

// calculate the correlation coefficient between two vector
double correlationCoefficient(std::vector<double> &X, std::vector<double> &Y);

// calculate the assortativity based on the pearson correlation coefficient
double getPearsonCorr(Graph& G);

// calculate the assortativity based on the spearman correlation coefficient
double getSpearmanCorr(Graph & G);

std::string tostr(int a);

std::string tostr(double a);

// write a graph to a file
void writeToFile(Graph &G, std::string file);


// write the results to file
void printResultsToFile(std::vector<std::vector<int> > degreeDists, std::vector<std::vector<double> > CdClusts, 
	std::vector<double> avgClusts, std::vector<double> globalClusts, std::vector<double> assorts, 
	std::vector<double> spearmans,std::string path, std::string modelname,std::vector<std::string> params);




#endif

