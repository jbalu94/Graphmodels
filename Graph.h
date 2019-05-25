#ifndef Graph_H
#define Graph_H

#include <vector>

struct Node
{
	int ind;
	std::vector<double> pos;
	Node() {}
	Node(int i, std::vector<double> p) { ind = i; pos = p; }
};

struct Edge
{
	Node u;
	Node v;
	Edge() {}
	Edge(Node u0, Node v0) { u = u0; v = v0; }

};

struct Graph
{
	std::vector<std::vector<Node> > edges;
	std::vector<Node> nodes;
	Graph() { edges.clear(); nodes.clear(); };
	void add_node(Node u)
	{
		edges.push_back(std::vector<Node>());
	}
	void add_edge(Node u, Node v)
	{
		edges[u.ind].push_back(v);
		edges[v.ind].push_back(u);
	}
	bool has_edge(Node u, Node v)
	{
		for (int i = 0; i < edges[u.ind].size(); ++i)
		{
			if (edges[u.ind][i].ind == v.ind) return true;
		}
		return false;
	}
};

#endif