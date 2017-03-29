#ifndef DHW_GEN_DEFORMATION_GRAPH_H_
#define DHW_GEN_DEFORMATION_GRAPH_H_

#include "fast_marching.h"
#include "deformation_graph.h"
#include <iostream>
#include <random>

double gen_DeformationGraph(DeformationGraph &dg,
	FastMarchingData &fmdata, const Eigen::MatrixXd &V,
	const Eigen::MatrixXi &F, const double radious_coef,
	const int node_nodes_num = 4, const int v_nodes_num = 6) 
{	
	vector<int> start_points(1);
	int nV = V.rows();

	FM_Float distmax = 0;

	Eigen::Vector3d m = V.colwise().minCoeff();
	Eigen::Vector3d M = V.colwise().maxCoeff();
	distmax = 3 * dg.nodes.size() / ((M(0) - m(0))*(M(1) - m(1))*(M(2) - m(2)));  //because a*b*c/x^3 = num £¬so x = sqrt3(...)

	if (fmdata.option.distmax > distmax)
		fmdata.option.distmax = distmax;

	//	fmdata.PerformFastMarching(V, F, dg.nodes); //do we need this to recompute from nodes?

	distmax = 0;
#ifdef AVARAGE_SEED_DISTANCE
	for (int i = 0; i < fmdata.seed_points.size(); i++) {
		int seed = fmdata.seed_points[i];
		FM_Float seed_max_dist = 0;
		for (int j = 0; j < V.rows(); j++) {
			if (fmdata.vertices[j].front == seed && fmdata.vertices[j].distance > seed_max_dist) {
				seed_max_dist = fmdata.vertices[j].distance;
			}
		}
		distmax += seed_max_dist;
	}
	distmax /= fmdata.seed_points.size();
#else
	for (int j = 0; j < V.rows(); j++) {
		if (fmdata.vertices[j].distance<FM_INFINITE && fmdata.vertices[j].distance > distmax) {
			distmax = fmdata.vertices[j].distance;
		}
	}
#endif

	std::cout << "the everage mad_dist distance from a seed is:" << distmax << "\n";

	//recomput the distance of a vertex  to a seed vertex and the distance no more than distmax
	fmdata.option.distmax = distmax*radious_coef;

	if (fmdata.option.bound) {
		delete[] fmdata.option.bound;
		fmdata.option.bound = 0;
	}
	dg.node_radious = fmdata.option.distmax;
	fmdata.seed_points.clear();

	dg.node_nodes.clear();
	dg.v_nodes.clear();
	dg.GNode_gs.resize(dg.nodes.size(), 3);  // node positions
	dg.node_nodes.resize(dg.nodes.size()); //node cnonectivity
	dg.v_nodes.resize(nV); // adj nodes and cooresponding weights of a vertex v 

	vector<int> seed_flag(nV, -1);
	int seed_v;
	for (int i = 0; i < dg.nodes.size(); i++) {
		seed_v = dg.nodes[i];
		seed_flag[seed_v] = i;  //seed_flag[seed_v] = i; the index in nodes of the seed vertex seed_v
		dg.GNode_gs.row(i) = V.row(seed_v);
	}
	for (int i = 0; i<dg.nodes.size(); i++) {
		int seed_v = dg.nodes[i];
		start_points[0] = seed_v;
		for (int j = 0; j < nV; j++) {
			fmdata.vertices[j].reset(); //	vertices[i].front = -1;
		}
		std::cout << "recompute distance from a seed vertex:" << start_points[0] << "\n";
		fmdata.PerformFastMarching(V, F, start_points);

		//costruct deformation graph
		std::cout << "costruct deformation graph:" << seed_v << "\n";
		for (int j = 0; j < nV; j++) {
			if (fmdata.vertices[j].distance < fmdata.option.distmax) {
				if (j != seed_v&&seed_flag[j] >= 0)
					dg.node_nodes[i].push_back(make_pair(seed_flag[j], fmdata.vertices[j].distance));
				//		node_nodes[i].push_back(seed_flag[j]);
				dg.v_nodes[j].push_back(make_pair(i, fmdata.vertices[j].distance));
			}
		}
	}
	if (dg.node_nodes.size()>node_nodes_num)
		reduce(dg.node_nodes, node_nodes_num);// 4);
	if (dg.v_nodes.size()>v_nodes_num)
		reduce(dg.v_nodes, v_nodes_num);// 6);
#if 1
	for (int j = 0; j < nV; j++)
		std::cout << dg.v_nodes[j].size() << "\n";
#endif
	return distmax;

}
double gen_DeformationGraph(DeformationGraph &dg, FastMarchingData &fmdata, const Eigen::MatrixXd &V,
	const Eigen::MatrixXi &F, const int num = 100, const double radious_coef = 2.1,
	const int node_nodes_num = 4, const int v_nodes_num = 6) 
{
	vector<int> start_points;
	std::random_device rd;
	std::mt19937 gen(rd());
	int nV = V.rows();
	std::uniform_int_distribution<> dis(0, nV - 1);
	for (int i = 0; i < 1; ++i)
		start_points.push_back(dis(gen));

	//	fmdata.option.iter_max =   10000;
	std::pair<int, FM_Float> p = fmdata.FarthestPointSampling(V, F, start_points, num);

	//recomput the distance of a vertex  to a seed vertex and the distance no more than distmax
	dg.nodes = fmdata.seed_points;   // node vertices
	return gen_DeformationGraph(dg,fmdata, V, F, radious_coef); // average_distance*radious_coef


}

double gen_DeformationGraph(DeformationGraph &dg, const Eigen::MatrixXd &V,
	const Eigen::MatrixXi &F, const int node_NUM = 100, const double radious_coef = 2.1,
	const int node_nodes_num = 4, const int v_nodes_num = 6)
{
	FastMarchingData fmdata;
	fmdata.option.iter_max = 1000;
	fmdata.PrepareFastMarching(V, F);	
	return gen_DeformationGraph(dg, fmdata, V, F, node_NUM, radious_coef, node_nodes_num, v_nodes_num);
	
}
#endif