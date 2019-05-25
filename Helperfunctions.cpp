
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <cmath>
#include "Graph.h"

#include "Helperfunctions.h"


double min(double a, double b) { return a < b ? a : b; }
double max(double a, double b) { return a < b ? b : a; }


//distance in the torus metric 
double distTorus(std::vector<double> a, std::vector<double>  b)
{
	double dist = 0;
	for (int i = 0; i < a.size(); ++i) dist += pow(min(min(fabs(a[i] - b[i]), fabs(1.0 - (a[i] - b[i]))), fabs(-1.0 - (a[i] - b[i]))), 2);
	return pow(dist, 0.5);
}


// dot product of two vector
double dot(std::vector<double> a, std::vector<double> b)
{
	double sum = 0;
	for (int i = 0; i < a.size(); ++i)
	{
		sum += a[i] * b[i];
	}
	return sum;
}


// calculate the angular between two vector
double getAng(std::vector<double> a, std::vector<double> b)
{
	return fabs(acos(dot(a, b) / (pow(dot(a, a), 0.5)*pow(dot(b, b), 0.5))));
}



// generate a random floating number between [0,1]
double fRand(double fMin, double fMax)
{
	double f = (double)rand() / RAND_MAX;
	return fMin + f * (fMax - fMin);
}



// rankify a vector
std::vector<double> rankify(std::vector<double> & X)
{

	int N = X.size();

	// Rank Vector 
	std::vector<double> Rank_X(N);

	for (int i = 0; i < N; i++)
	{
		int r = 1, s = 1;

		// Count no of smaller elements 
		// in 0 to i-1 
		for (int j = 0; j < i; j++) {
			if (X[j] < X[i]) r++;
			if (X[j] == X[i]) s++;
		}

		// Count no of smaller elements 
		// in i+1 to N-1 
		for (int j = i + 1; j < N; j++) {
			if (X[j] < X[i]) r++;
			if (X[j] == X[i]) s++;
		}

		// Use Fractional Rank formula 
		// fractional_rank = r + (n-1)/2 
		Rank_X[i] = r + (s - 1) * 0.5;
	}

	// Return Rank Vector 
	return Rank_X;
}



// generating a random point on a sphere
std::vector<double> generateRandomToSphere(double r)
{
	// assuming that radius = 1/(2*sqrt(pi))
	double x = fRand(-r, r);
	double y = fRand(-r, r);
	double z = fRand(-r, r);
	double l = pow(x*x + y * y + z * z, 0.5);
	if (l < r)
	{
		double scalar = r / l;
		return std::vector<double>{x*scalar, y*scalar, z*scalar};
	}
	else return generateRandomToSphere(r);
}


// calculate degree distribution
std::vector<int> getDegreeDist(Graph &G)
{
	int n = G.nodes.size();
	std::vector<int> counters(G.nodes.size());
	for (int i = 0; i < n; i++)
	{
		if (G.edges[i].size() < counters.size()) counters[G.edges[i].size()]++;
	}
	return counters;
}


// calculate global clustring coefficient
double getGlobalClustering(Graph &G)
{

	int cherries = 0;
	int triangles = 0;
	std::set<std::vector<int> > tris;
	for (int i = 0; i < G.nodes.size(); ++i)
	{
		std::set<int> neighs;
		for (int j1 = 0; j1 < G.edges[i].size(); ++j1)
		{
			if (neighs.find(G.edges[i][j1].ind) == neighs.end() && G.edges[i][j1].ind != G.nodes[i].ind) neighs.insert(G.edges[i][j1].ind);
			for (int j2 = 0; j2 < G.edges[i].size(); ++j2)
			{
				Node u = G.edges[i][j1];
				Node v = G.edges[i][j2];
				if (u.ind < v.ind && u.ind > G.nodes[i].ind && v.ind > G.nodes[i].ind && G.has_edge(u, v) &&
					tris.find({ G.nodes[i].ind ,u.ind,v.ind }) == tris.end())
				{
					triangles++;
					tris.insert({ G.nodes[i].ind ,u.ind,v.ind });
				}
			}
		}
		cherries += neighs.size()*(neighs.size() - 1) / 2;
	}
	//cout << endl << "triangles: " << triangles << endl;
	//cout << "cherries: " << cherries << endl;
	std::cout << std::endl << std::endl << "global clustering coeff: " << double(3 * triangles / double((cherries))) << std::endl;

	return double(3 * triangles / double((cherries)));
}


// calculate avg clustering coefficient
double getAvgClustering(Graph &G)
{
	double avgCoeff = 0;
	int notCounted = 0;
	for (int i = 0; i < G.nodes.size(); ++i)
	{
		int triangles = 0;
		for (int j1 = 0; j1 < G.edges[i].size(); ++j1)
		{
			for (int j2 = j1 + 1; j2 < G.edges[i].size(); ++j2)
			{
				Node u = G.edges[i][j1];
				Node v = G.edges[i][j2];
				if (u.ind != v.ind && u.ind != G.nodes[i].ind && v.ind != G.nodes[i].ind &&G.has_edge(u, v)) triangles++;
			}
		}
		if (G.edges[i].size() > 1)
		{
			avgCoeff += double(triangles / double((G.edges[i].size()*(G.edges[i].size() - 1) / 2)));
		}
		else
		{
			notCounted++;
		}
	}
	avgCoeff /= double(G.nodes.size() - notCounted);
	std::cout << "Average clustering coefficient:" << avgCoeff << std::endl;

	return avgCoeff;
}


// calculate C(d) function, which is the mean of the local clustering coefficient which has 
std::vector<double> getCd(Graph& G)
{
	std::map<int, std::vector<double> > C;
	int n = G.nodes.size();

	double avgCoeffs = 0;
	int notCounted = 0;
	for (int i = 0; i < G.nodes.size(); ++i)
	{
		int triangles = 0;
		for (int j1 = 0; j1 < G.edges[i].size(); ++j1)
		{
			for (int j2 = j1 + 1; j2 < G.edges[i].size(); ++j2)
			{
				Node u = G.edges[i][j1];
				Node v = G.edges[i][j2];
				if (u.ind != v.ind && u.ind != G.nodes[i].ind && v.ind != G.nodes[i].ind &&G.has_edge(u, v)) triangles++;
			}
		}
		if (G.edges[i].size() > 1)
		{
			int deg = G.edges[i].size();
			double coeff = double(triangles / double((G.edges[i].size()*(G.edges[i].size() - 1) / 2)));
			C[deg].push_back(coeff);
			//avgCoeff += double(triangles / double((G.edges[i].size()*(G.edges[i].size() - 1) / 2)));
		}
		else
		{
			notCounted++;
		}
	}
	std::vector<double> C_d(n, 0);
	for (int i = 0; i < n; ++i)
	{
		if (C[i].size() > 0)
		{
			double sum = 0;
			for (int j = 0; j < C[i].size(); ++j)
			{
				sum += C[i][j];
			}
			C_d[i] = double(sum) / double(C[i].size());
		}
	}
	return C_d;
}


// calculate the correlation coefficient between two vector
double correlationCoefficient(std::vector<double> &X, std::vector<double> &Y)
{
	int n = X.size();
	float sum_X = 0, sum_Y = 0,
		sum_XY = 0;
	float squareSum_X = 0,
		squareSum_Y = 0;

	for (int i = 0; i < n; i++)
	{
		// sum of elements of array X. 
		sum_X = sum_X + X[i];

		// sum of elements of array Y. 
		sum_Y = sum_Y + Y[i];

		// sum of X[i] * Y[i]. 
		sum_XY = sum_XY + X[i] * Y[i];

		// sum of square of array elements. 
		squareSum_X = squareSum_X +
			X[i] * X[i];
		squareSum_Y = squareSum_Y +
			Y[i] * Y[i];
	}

	// use formula for calculating 
	// correlation coefficient. 
	double corr = (double)(n * sum_XY -
		sum_X * sum_Y) /
		sqrt((n * squareSum_X -
			sum_X * sum_X) *
			(n * squareSum_Y -
				sum_Y * sum_Y));

	return corr;
}

double getPearsonCorr(Graph& G)
{
	std::vector<double> x;
	std::vector<double> y;
	for (int i = 0; i < G.edges.size(); ++i)
	{
		for (int j = 0; j < G.edges[i].size(); ++j)
		{
			int u = i;
			int v = G.edges[i][j].ind;
			int deg1 = G.edges[u].size();
			int deg2 = G.edges[v].size();
			x.push_back(deg1);
			y.push_back(deg2);
		}
	}
	double corr = correlationCoefficient(x, y);
	std::cout << "Assortativity with pearson correlation coefficient: " << corr << std::endl;

	return corr;
}

double getSpearmanCorr(Graph & G)
{
	std::vector<double> x;
	std::vector<double> y;
	for (int i = 0; i < G.edges.size(); ++i)
	{
		for (int j = 0; j < G.edges[i].size(); ++j)
		{
			int u = i;
			int v = G.edges[i][j].ind;
			int deg1 = G.edges[u].size();
			int deg2 = G.edges[v].size();
			x.push_back(deg1);
			y.push_back(deg2);
		}
	}
	std::vector<double> ranked_x = rankify(x);
	std::vector<double> ranked_y = rankify(y);
	double corr = correlationCoefficient(ranked_x, ranked_y);
	std::cout << "Assortativity with spearman correlation coefficient: " << corr << std::endl;
	return corr;
}


std::string tostr(int a)
{
	std::stringstream ss;
	ss << a;
	return ss.str();
}
std::string tostr(double a)
{
	std::stringstream ss;
	ss << a;
	return ss.str();
}