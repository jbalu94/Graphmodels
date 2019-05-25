#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <random>
#include <fstream>


#include "Graph.h"
#include "Helperfunctions.h"
#include "HPAgraph.h"

#define PI 3.1415926



// calculate the hiperbolic distance between two points
double getHiperDist(std::vector<double> a, std::vector<double> b, double beta)
{
	double s = a[0];
	double t = b[0];
	double s_ang = a[1];
	double t_ang = b[1];

	double ang_dist = min(fabs(s_ang - t_ang), min(fabs(s_ang - 2 * PI - t_ang), fabs(s_ang + 2 * PI - t_ang)));

	double r_t = 2.0 * log(t);
	double r_s = 2.0 * log(s);
	double r_st = beta * r_s + (1 - beta)*r_t;

	return r_st + r_t + 2.0* log(ang_dist / 2.0);
}


// generate the ang of the new node
double generateHPAang(Graph &G, int t, double beta, double delta, std::mt19937& gen)
{
	std::vector<double> randAngs(t);

	double r_t = 2 * log(t);

	for (int i = 0; i < t; ++i)
	{
		randAngs[i] = fRand(0, 2 * PI);
	}

	std::vector<double> closePointNum(t, delta);
	for (int i = 0; i < t; ++i)
	{
		for (int j = 0; j < t - 1; ++j)
		{
			if (getHiperDist(std::vector<double>{r_t, randAngs[i]}, G.nodes[j].pos, beta) < r_t)
			{
				closePointNum[i] += 1;
			}
		}
	}

	std::discrete_distribution<int> dist(std::begin(closePointNum), std::end(closePointNum));
	int rand = dist(gen);
	return randAngs[rand];
}

Graph generateHPA(int n, int m, double beta, double delta, std::mt19937& gen)
{

	Graph G = Graph();
	G.nodes.resize(n);
	G.edges.resize(n);

	double firstNodeAng = fRand(0, 2 * PI);
	Node firstnode = Node(0, std::vector<double>{1.0, firstNodeAng});
	G.nodes[0] = firstnode;

	for (int i = 1; i < n; ++i)
	{
		double ang = generateHPAang(G, i, beta, delta, gen);
		G.nodes[i] = Node(i, std::vector<double>{double(i + 1), ang});
		if (i <= m)
		{
			for (int j = 0; j < i; ++j)
			{
				G.add_edge(G.nodes[i], G.nodes[j]);
			}
		}
		else
		{
			std::vector<std::pair<double, int> > distsFromCurrentNode(i);
			for (int j = 0; j < i; ++j)
			{
				double dist = getHiperDist(G.nodes[j].pos, G.nodes[i].pos, beta);
				distsFromCurrentNode[j] = std::make_pair(dist, j);
			}
			std::sort(distsFromCurrentNode.begin(), distsFromCurrentNode.end());
			for (int j = 0; j < m; ++j)
			{
				G.add_edge(G.nodes[i], G.nodes[distsFromCurrentNode[j].second]);
			}
		}
	}
	return G;
}

// simulate HPA graphs with certain parameters, and save it to files
void simulateHPAandSave(std::string path, int T, std::vector<int> N, std::vector<int> M,
	std::vector<double> betas, std::vector<double> deltas, std::mt19937& gen)
{
	int I = 0;
	for (auto n : N)
	{
		for (auto m : M)
		{
			for (auto beta : betas)
			{
				
				for (auto delta : deltas)
				{
					std::vector<std::vector<int> > degreeDists(T);
					std::vector<double> avgClusts(T);
					std::vector<double> globalClusts(T);
					std::vector<std::vector<double> > CdClusts(T);
					std::vector<double> assorts(T);
					std::vector<double> spearmans(T);

					for (int t = 0; t < T; ++t)
					{
						std::cout << std::endl << std::endl << ++I << ":   n=" << n << "  m=" << m << "   beta="
							<< beta << "   delta=" << delta << "   t=" << t;

						Graph G = generateHPA(n, m, beta, delta, gen);
						std::vector<int> degreeDist = getDegreeDist(G);
						std::vector<double> cd = getCd(G);
						double globalClust = getGlobalClustering(G);
						double dn = double(n);

						double avgClust = getAvgClustering(G);
						double pearson = getPearsonCorr(G);
						double spearman = getSpearmanCorr(G);


						degreeDists[t] = degreeDist;
						CdClusts[t] = cd;
						avgClusts[t] = avgClust;
						globalClusts[t] = globalClust;
						assorts[t] = pearson;
						spearmans[t] = spearman;
					}

					std::string n_str = tostr(n);
					std::string m_str = tostr(m);
					std::string beta_str = tostr(beta);
					std::string delta_str = tostr(delta);

					std::vector<std::string>  params{ n_str,m_str,beta_str,delta_str };
					printResultsToFile(degreeDists, CdClusts, avgClusts, globalClusts, assorts, spearmans, path, "HPA", params);
				}
			
			}
		}
	}
}