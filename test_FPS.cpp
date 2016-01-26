#include <fstream>
#include <string>
#include <iostream>
#include <Eigen/Dense>
#include "fast_marching.h"
#include "deformation_graph.h"
#include <random>

using std::ifstream;
using std::ofstream;
using std::string;
using std::cout;

template <typename DerivedV, typename DerivedF>
inline int readOFF(
	const char * str,
	Eigen::PlainObjectBase<DerivedV>& V,
	Eigen::PlainObjectBase<DerivedF>& F)
{
	ifstream iFile(str);
	if (!iFile) return 1;
	string fT;
	char ch;
	int nV, nF, nN;
	iFile >> fT;

	if (fT.compare("OFF") != 0 && fT.compare("off") != 0) return 2;
	iFile >> nV >> nF >> nN;

	V.setConstant(nV, 3, 0);
	for (int i = 0; i < nV; i++)
		iFile >> V(i, 0) >> V(i, 1) >> V(i, 2);

	F.setConstant(nF, 3, 0);
	for (int i = 0; i < nF; i++)
		iFile >> ch >> F(i, 0) >> F(i, 1) >> F(i, 2);

#if 0

	for (int i = 0; i < nV; i++) {
		std::cout << V(i, 0) << " " << V(i, 1) << " " << V(i, 2) << std::endl;
	}
	for (int i = 0; i < nF; i++) {
		std::cout << F(i, 0) << " " << F(i, 1) << " " << F(i, 2) << std::endl;
	}

#endif
	return 0;
}

Eigen::MatrixXd V;
Eigen::MatrixXi F;
FastMarchingData fmdata;
bool saveFSP_Result(const char *file, FastMarchingData &fmdata) {
	ofstream oFile(file); 
	if (!oFile) return false;
	int nS = fmdata.seed_points.size();	
	int nV = fmdata.vertices.size();
	oFile << nV<<" "<<nS<<" "<< endl;
	for (auto v : fmdata.seed_points)
		oFile << v <<endl;
	for (auto v : fmdata.vertices) {
		oFile << v.front<<" " << v.state << " " <<v.distance << endl;
	}
	oFile.close();
	return true;
}
int test_FPS(int argc, char *argv[]) {
	if (argc < 2) {
		std::cout << "lack of mesh file and result file!\n";
		std::cout << "run this program in format: test_FPS meshfile resultfile \n";
		return 1;
	}
	if (readOFF(argv[1], V, F) != 0) {
		std::cout << "read the file " << argv[1] << " failed";
		return 2;
	}
	vector<int> start_points;
	fmdata.option.iter_max = 1000;

	std::cout << "PrepareFastMarching...\n";
	fmdata.PrepareFastMarching(V, F);

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<> dis(0, V.rows() - 1);
	int n = 1;
	for (int i = 0; i < n; ++i)
		start_points.push_back(dis(gen));
	start_points[0] = 18842;
	//  start_points.push_back(18842);
	// start_points.push_back(9417);
	//	fmdata.PerformFastMarching(V, F, start_points);
	std::cout << "start FarthestPointSampling...\n";
	fmdata.FarthestPointSampling(V, F, start_points, 20);

	if (argc >= 3) {
		std::cout << "save to file "<< argv[2]<< "....\n";
		saveFSP_Result(argv[2], fmdata);
	}
	return 0;
}

int main(int argc, char *argv[]) {
	return test_FPS(argc,argv);
}