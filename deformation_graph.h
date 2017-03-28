
#ifndef DHW_DEFORMATION_GRAPH_H
#define DHW_DEFORMATION_GRAPH_H

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>
#include <fstream>
#include <algorithm>
#include <functional>
#include "fast_marching.h"
using namespace std;

struct DeformationGraph {
	/*
	struct Graph_Node {
		//	double g[3];
		//	double R[9];
		//	double t[3];
		vector<int> AdjNodes;
	};*/
	Eigen::MatrixXd GNode_gs;
	vector<int> nodes;
	//vector<Graph_Node> GNodes;
	vector< vector<pair <int, double>> > node_nodes;// node 邻接的GraphNode及权系数
	vector< vector<pair <int, double>> >  v_nodes;  //顶点v 邻接的GraphNode及权系数
//	vector< vector<double> > m_wi_vj;//vj 邻接的GraphNode的权系数
//	vector< vector<int> > m_gi_vj; //vj 邻接的GraphNode.
	double global_g[3];

	double node_radious;

	bool save_DeformationGraph(const std::string filename);
	bool read_DeformationGraph(const std::string filename);
};

inline void reduce(vector< vector<pair <int, double>> > &adj_nodes, const int maxNum) {
	for (int i = 0; i < adj_nodes.size(); i++)
		if (adj_nodes[i].size()>maxNum) {
			std::sort(adj_nodes[i].begin(), adj_nodes[i].end(), 
			      [](const pair <int, double> &a, const pair <int, double> &b) {
					  return a.second < b.second;
				  });
			adj_nodes[i].resize(maxNum);
		}
}

#endif // !define 