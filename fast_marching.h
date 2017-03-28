#ifndef DHW_FAST_MARCHING_H
#define DHW_FAST_MARCHING_H

#include "igl/igl_inline.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>
#include <memory>
using namespace std;
namespace igl
{

}
typedef double FM_Float;
const FM_Float FM_INFINITE = 1e15;


struct FastMarchingOption {
	vector<int> start_points = vector<int>();
	vector<double> start_values = vector<double>();
	vector<int> _end_points = vector<int>();
	FM_Float distmax = 1e9;
	int iter_max = 100000000;//-1,
	FM_Float* weight = 0;	// weight
	FM_Float* heuristic = 0;	// heuristic
	FM_Float* bound = 0;	// bound on current distance
	bool bUseUnfolding_ = true;
};

class FastMarchingVertex {
public:
	enum FastMarchingVertexState
	{
		kFar,
		kAlive,
		kDead
	};
	int vid;
	int front=-1;
	FM_Float distance = FM_INFINITE;
	FastMarchingVertexState state = kFar;
	bool stopped = false;

	struct Front_Distance {
		int i; FM_Float value= FM_INFINITE;
		Front_Distance(int ii, FM_Float v) :i(ii), value(v) {}
	};
	struct FrontOverlapInfo {
		vector<Front_Distance> fronts;
	};
	std::shared_ptr<FrontOverlapInfo> overlapInfo;

	void reset() {
		front = -1; distance = FM_INFINITE; state = kFar; stopped = false; overlapInfo.reset();	}
};

typedef FastMarchingVertex::FastMarchingVertexState VertexState;

typedef std::vector<class FastMarchingVertex*> VertexPtrVector;



struct FastMarchingData {
	// mesh connectivity
	//std::vector<bool> V_border; //Eigen::MatrixXi E, E2F, F2E;
	std::vector<std::vector<int> > VF, VFi, VV;
	Eigen::MatrixXi TT, TTi;
	
	//Fast Matching data
	std::vector<struct FastMarchingVertex> vertices;

	/*
	// mesh V,F
	Eigen::MatrixXd &V;
	Eigen::MatrixXi &F;
	FastMarchingData(Eigen::MatrixXd &mV, Eigen::MatrixXi &mF) :V(mV), F(mF) {};
	*/

	FastMarchingOption option;

	//result
	vector<double> distance;	// distance
	vector<double> state;	// state
	vector<double> nearest_neighbor;

	vector<int> seed_points;
	bool BuildConnectivity(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F);
	bool PerformFastMarching(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
		const vector<int> &start_points);
	bool PrepareFastMarching(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F);
	std::pair<int, FM_Float> FarthestPointSamplingStep(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F);
	std::pair<int, FM_Float> FarthestPointSampling(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
		const vector<int> &start_points,int num = 0);
	
private:
	int nbr_iter = 0;
	void NewDeadVertexCallback_(FastMarchingVertex* pV=0) { }
	bool VertexInsersionCallback_(FastMarchingVertex* pV, double new_dist) { 
		if (nbr_iter > option.iter_max) return false; 
		if (option.bound&&new_dist > option.bound[pV->vid]) return false;
		return true;
	}
	bool frontOverlap(FastMarchingVertex* pNewVert,int new_front,FM_Float new_distance);
	bool toStop(FastMarchingVertex* pVert);
};

#ifndef IGL_STATIC_LIBRARY
//#include "fast_marching.cpp"
#endif

#endif
