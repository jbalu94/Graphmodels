
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <random>
#include <fstream>


#include "GPAgraph.h"
#include "Graph.h"
#include "Helperfunctions.h"

#define PI 3.1415926

double getMinBucketDistEstimate(std::vector<std::vector<std::vector<std::vector<double> > > > &bucketSarkok,
	std::vector<double> p, int i, int j, double unit)
{

	std::vector<double> D(4);

	D[0] = distTorus(p, bucketSarkok[i][j][0]);
	D[1] = distTorus(p, bucketSarkok[i][j][1]);
	D[2] = distTorus(p, bucketSarkok[i][j][2]);
	D[3] = distTorus(p, bucketSarkok[i][j][3]);

	sort(D.begin(), D.end());
	if (D[0] * D[0] + unit * unit <= D[1] * D[1]) return D[0];
	return D[0] * 0.708;

}


double getMaxBucketDist(std::vector<std::vector<std::vector<std::vector<double> > > > &bucketSarkok, std::vector<double> p, int i, int j)
{
	double d1 = distTorus(p, bucketSarkok[i][j][0]);
	double d2 = distTorus(p, bucketSarkok[i][j][1]);
	double d3 = distTorus(p, bucketSarkok[i][j][2]);
	double d4 = distTorus(p, bucketSarkok[i][j][3]);
	return max(d1, max(d2, max(d3, d4)));
}


// It calculates the expected number of nodes with degree k after inserting n nodes to the PA model 
double getProb(int n, int m, double delta, int k)
{
	if (delta == 0)
	{
		return double(n * 2 * m*(m + 1) / (k*(k + 1)*(k + 2)));
	}

	double md = double(m);
	return n * (2.0 + delta / md)*(tgamma(k + delta)*tgamma(md + 2 + delta + delta / md)) / 
		(tgamma(md + delta)*tgamma(k + 3 + delta + delta / md));
}



// calculates the value of F function, in the case when F = F0
double getFvalue0(double r, double u)
{
	return 1;
}


// calculates the value of F function, in the case when F = F1
double getFvalue1(double r, double u)
{
	if (fabs(u) < r) return 1.0;
	else return 0;
}


// calculates the value of F function, in the case when F = F2
double getFvalue2(double beta, double u, int n, double theta = 0.25)
{
	return pow(max(u, pow(double(n), -theta)), -beta);
}


// calculates the value of F function, in the case when F = F3
double getFvalue3(double beta, double u)
{
	return exp(-beta * u);
}



// calculate the integral which appered in the GPA model, and was denoted by I
double getI(double beta, int n)
{
	std::vector<double> x(10000);
	std::vector<double> y(10000);
	x[0] = 0.00001;
	double bound = pow(double(n), -0.25);
	for (int i = 1; i < x.size(); ++i) x[i] = x[i - 1] + PI / double(x.size());
	for (int i = 0; i < y.size(); ++i)
	{
		y[i] = bound < x[i] ? pow(x[i], -beta)*sin(x[i]) : pow(bound, -beta)*sin(bound);
	}

	double int_d = 0.0;
	double int_u = 0.0;
	for (int i = 0; i < x.size() - 1; ++i)
	{
		int_d += y[i] * (x[i + 1] - x[i]);
		int_u += y[i + 1] * (x[i + 1] - x[i]);
	}
	return int_d;
}



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

Graph generateGPA1(int n, int m, double delta, double r, int k, std::mt19937& gen)
{
	bool toosmall = false;
	if (double(1.0 / k)*0.707 > r)
	{
		std::cout << "WARNING!! it can be possible that inside the bucket but not good" << std::endl;
		toosmall = true;
	}

	std::vector<std::vector<std::vector<int> > > buckets(k, std::vector<std::vector<int> >(k, std::vector<int>()));
	std::vector<std::vector<std::vector<std::vector<double> > > > bucketSarkok
	(k, std::vector<std::vector<std::vector<double> > >(k, std::vector<std::vector<double> >(4)));

	double unit = 1.0 / k;
	for (int i = 0; i < k; i++)
	{
		for (int j = 0; j < k; ++j)
		{
			bucketSarkok[i][j][0] = std::vector<double>{ i*unit, j*unit };
			bucketSarkok[i][j][1] = std::vector<double>{ i*unit, (j + 1)*unit };
			bucketSarkok[i][j][2] = std::vector<double>{ (i + 1)*unit, j*unit };
			bucketSarkok[i][j][3] = std::vector<double>{ (i + 1)*unit, (j + 1)*unit };
		}
	}

	Graph G = Graph();

	G.nodes.resize(n);
	G.edges.resize(n);

	for (int i = 0; i < n; ++i)
	{
		double X = ((double)rand() / (double)RAND_MAX);
		double Y = ((double)rand() / (double)RAND_MAX);
		G.nodes[i] = Node(i, std::vector<double>{X, Y});
	}

	int u = int(G.nodes[0].pos[0] / unit);
	int v = int(G.nodes[0].pos[1] / unit);
	if (u == k) u--;
	if (v == k) v--;
	buckets[u][v].push_back(0);
	G.edges[0].push_back(Node(0, std::vector<double>{ 0, 0 }));

	for (int i = 1; i < n; ++i)
	{
		int cc = 0;
		//if (i % 1000 == 0) cout << i << " ";
		std::vector<double> actPos = G.nodes[i].pos;
		std::vector<int> closePoints;
		for (int i1 = 0; i1 < k; ++i1)
		{
			for (int i2 = 0; i2 < k; i2++)
			{
				double mindist = getMinBucketDistEstimate(bucketSarkok, actPos, i1, i2, unit);
				double maxdist = getMaxBucketDist(bucketSarkok, actPos, i1, i2);

				if (maxdist <= r && !toosmall)
				{
					closePoints.insert(closePoints.end(), buckets[i1][i2].begin(), buckets[i1][i2].end());
				}
				else if (mindist <= r || toosmall)
				{
					cc++;
					for (int i3 = 0; i3 < buckets[i1][i2].size(); ++i3)
					{
						if (distTorus(actPos, G.nodes[buckets[i1][i2][i3]].pos) <= r)
						{
							closePoints.push_back(buckets[i1][i2][i3]);
						}
					}
				}
			}
		}

		std::vector<double> sizes(closePoints.size());
		for (int j = 0; j < closePoints.size(); ++j)
		{
			sizes[j] = G.edges[closePoints[j]].size() + delta;
		}
		closePoints.push_back(i);
		sizes.push_back(1 + delta);
		for (int z = 0; z < m; ++z)
		{
			std::discrete_distribution<int> dist(std::begin(sizes), std::end(sizes));
			int rand = dist(gen);
			G.add_edge(G.nodes[i], G.nodes[closePoints[rand]]);
			if (rand != closePoints.size() - 1) sizes[rand] += 1;
		}

		int u = int(G.nodes[i].pos[0] / unit);
		int v = int(G.nodes[i].pos[1] / unit);
		if (u == k) u--;
		if (v == k) v--;
		buckets[u][v].push_back(i);
	}

	return G;
}





// generating graph from the model when F = F2
// this is on the sphere, with the angular metric
Graph generateGPA2(int n, int m, double alpha, double delta, double beta, std::mt19937& gen)
{
	// first, generate the points
	Graph G = Graph();
	G.nodes.resize(n);
	G.edges.resize(n);
	for (int i = 0; i < n; ++i)
	{
		G.nodes[i].ind = i;
		G.nodes[i].pos = generateRandomToSphere(0.5 / sqrt(PI));
	}

	// generate edges
	for (int i = 1; i <= max(1, m / 2); ++i)
	{
		G.add_edge(G.nodes[0], G.nodes[0]);
	}


	for (int i = 1; i < n; ++i)
	{
		std::vector<double> weights(i + 1, 0);
		double expected = getI(beta, n)*(2 * m + delta)*i;
		double summa = 0;
		for (int j = 0; j < i; ++j)
		{
			weights[j] = (G.edges[j].size() + delta)*getFvalue2(beta, getAng(G.nodes[j].pos, G.nodes[i].pos), n);
			summa += weights[j];
		}
		if (summa == 0) weights[i] = 1;
		else weights[i] = summa < alpha*expected / 2.0 ? (1 - summa / (alpha*expected / 2.0))*summa : 0;

		for (int z = 0; z < m; ++z)
		{
			std::discrete_distribution<int> dist(std::begin(weights), std::end(weights));
			int rand = dist(gen);
			G.add_edge(G.nodes[i], G.nodes[rand]);
		}
	}
	return G;
}


void simulateGPA2andSave(std::string path, int  T, std::vector<int> N, std::vector<int> M, std::vector<double> betas,
	std::vector<double> alphas, std::vector<double> deltas, std::mt19937& gen)
{
	int I = 0;
	for (auto n : N)
	{
		for (auto m : M)
		{
			for (auto beta : betas)
			{
				for (auto alpha : alphas)
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
							std::cout << ++I << ":   n=" << n << "  m=" << m << "   beta="
								<< beta << "  alpha=" << alpha << "   delta=" << delta << "   t=" << t << std::endl;

							Graph G = generateGPA2(n, m, alpha, delta, beta, gen);
							std::vector<int> degreeDist = getDegreeDist(G);
							std::vector<double> cd = getCd(G);
							double globalClust = getGlobalClustering(G);
							double dn = double(n);
							std::cout << "Expected triangles: " << (m - 1)*m*(m + 1)*log(dn)*log(dn)*log(dn) / 48.0 << std::endl;
							std::cout << "Expected cherries: " << m * (m + 1)*n*log(n) / 2.0 << std::endl;
							std::cout << "Expected global clustering coefficient: " << (m - 1) / 8.0*log(dn)*log(dn) / dn << std::endl;
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
						std::string alpha_str = tostr(alpha);

						std::string delta_str = tostr(delta);

						std::string degreeDistFile = path + "/degreeDist_" + n_str + "_" + m_str + "_" + beta_str + "_" + alpha_str + "_" + delta_str + ".txt";
						std::string avgClustFile = path + "/avgClust_" + n_str + "_" + m_str + "_" + beta_str + "_" + alpha_str + "_" + delta_str + ".txt";
						std::string globalClustFile = path + "/globalClust_" + n_str + "_" + m_str + "_" + beta_str + "_" + alpha_str + "_" + delta_str + ".txt";
						std::string CdClustFile = path + "/CdClust_" + n_str + "_" + m_str + "_" + beta_str + "_" + alpha_str + "_" + delta_str + ".txt";
						std::string assortFile = path + "/assort_" + n_str + "_" + m_str + "_" + beta_str + "_" + alpha_str + "_" + delta_str + ".txt";
						std::string spearmanFile = path + "/spearman_" + n_str + "_" + m_str + "_" + beta_str + "_" + alpha_str + "_" + delta_str + ".txt";

						std::ofstream f1(degreeDistFile.c_str());
						for (int i = 0; i < T; ++i)
						{
							for (int j = 0; j < n; ++j)
							{
								f1 << degreeDists[i][j] << " ";
							}
							f1 << std::endl;
						}
						f1.close();

						std::ofstream f2(avgClustFile.c_str());
						for (int i = 0; i < T; ++i)
						{
							f2 << avgClusts[i] << " ";
						}
						f2.close();

						std::ofstream f3(globalClustFile.c_str());
						for (int i = 0; i < T; ++i)
						{
							f3 << globalClusts[i] << " ";
						}
						f3.close();

						std::ofstream f4(assortFile.c_str());
						for (int i = 0; i < T; ++i)
						{
							f4 << assorts[i] << " ";
						}
						f4.close();

						std::ofstream f5(CdClustFile.c_str());
						for (int i = 0; i < T; ++i)
						{
							for (int j = 0; j < n; ++j)
							{
								f5 << CdClusts[i][j] << " ";
							}
							f5 << std::endl;
						}
						f5.close();

						std::ofstream f6(spearmanFile.c_str());
						for (int i = 0; i < T; ++i)
						{
							f6 << spearmans[i] << " ";
						}
						f6.close();
					}
				}
			}
		}
	}
}

