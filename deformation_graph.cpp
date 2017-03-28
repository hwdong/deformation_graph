
#include "deformation_graph.h"
#include <iostream>

#include <random>

bool DeformationGraph::save_DeformationGraph(const std::string filename) {
	ofstream oFile(filename);
	if (!oFile) return false;
	const int gn = GNode_gs.rows();
	oFile << v_nodes.size()<<" "<< gn <<" "<< node_radious << endl;
	for (int i = 0; i<gn; i++)
		oFile << nodes[i]<<" "<<GNode_gs(i, 0) << " " << GNode_gs(i, 1) << " " << GNode_gs(i, 2) << endl;

	int J;
	for (int i = 0; i<node_nodes.size(); i++) {
		oFile << node_nodes[i].size() << " ";
		for (int j = 0; j<node_nodes[i].size(); j++) {			
			oFile << node_nodes[i][j].first<<" " << node_nodes[i][j].second << " ";
		}
		oFile << "\n";
	}
	for (int i = 0; i < v_nodes.size(); i++) {
		oFile << v_nodes[i].size() << " ";
		for (int j = 0; j<v_nodes[i].size(); j++) {
			oFile << v_nodes[i][j].first << " " << v_nodes[i][j].second << " ";
		}
		oFile << "\n";
	}
	oFile.close();
	std::cout << "deformation graph has been saved to file" << filename << "\n";
	return true;
}

bool DeformationGraph::read_DeformationGraph(const std::string filename) {
	ifstream iFile(filename);
	if (!iFile) return false;
	int vn,gn;
	iFile >>vn>> gn>> node_radious;
	nodes.resize(gn);
	GNode_gs.resize(gn,3);
	node_nodes.resize(gn);
	v_nodes.resize(vn);
	for (int i = 0; i<gn; i++)
		iFile >> nodes[i] >> GNode_gs(i, 0)>> GNode_gs(i, 1)>> GNode_gs(i, 2) ;

	int n;
	for (int i = 0; i<gn; i++) {
		iFile >> n; node_nodes[i].resize(n);
		for (int j = 0; j<n; j++) {
			iFile >> node_nodes[i][j].first >> node_nodes[i][j].second;
		//	iFile >> node_nodes[i][j];
		}
	}
	for (int i = 0; i < vn; i++) {
		iFile >> n;  v_nodes[i].resize(n) ;
		for (int j = 0; j<n; j++) {
			iFile >> v_nodes[i][j].first >> v_nodes[i][j].second ;
		}
	}
	iFile.close();
	return true;
}


