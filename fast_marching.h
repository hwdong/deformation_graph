
#ifndef DHW_FAST_MARCHING_H
#define DHW_FAST_MARCHING_H

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>
#include <memory>
using namespace std;

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
	std::shared_ptr<FrontOverlapInfo> overlaoInfo;

	void reset() {
		front = -1; distance = FM_INFINITE; state = kFar; stopped = false; overlaoInfo.reset();	}
};

typedef FastMarchingVertex::FastMarchingVertexState VertexState;

typedef std::vector<FastMarchingVertex*> VertexPtrVector;



struct FastMarchingData {
	// mesh connectivity
	//std::vector<bool> V_border; //Eigen::MatrixXi E, E2F, F2E;
	std::vector<std::vector<int> > VF, VFi, VV;
	Eigen::MatrixXi TT, TTi;
	
	//Fast Matching data
	std::vector<FastMarchingVertex> vertices;

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
	int FarthestPointSamplingStep(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F);
	bool FarthestPointSampling(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
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

template <typename Index, typename IndexVector>
void vertex_vertex_adjacency(const Eigen::PlainObjectBase<Index>  & F,
	std::vector<std::vector<IndexVector> >& VV){
	VV.clear();
	VV.resize(F.maxCoeff() + 1);

	// Loop over faces
	for (int i = 0; i<F.rows(); i++){
		// Loop over this face
		for (int j = 0; j<F.cols(); j++){
			// Get indices of edge: s --> d
			int s = F(i, j);
			int d = F(i, (j + 1) % F.cols());
			VV.at(s).push_back(d);
			VV.at(d).push_back(s);
		}
	}
	// Remove duplicates
	for (int i = 0; i<(int)VV.size(); ++i){
		std::sort(VV[i].begin(), VV[i].end());
		VV[i].erase(std::unique(VV[i].begin(), VV[i].end()), VV[i].end());
	}
	
}

template <typename DerivedF, typename VFType>
void vertex_triangle_adjacency(	
	const Eigen::PlainObjectBase<DerivedF>& F,
	std::vector<std::vector<VFType> >& VF)
{
	int n = F.maxCoeff() + 1;
	VF.clear();	VF.resize(n);

	typedef typename DerivedF::Index Index;
	for (Index fi = 0; fi<F.rows(); ++fi)	{
		for (Index i = 0; i < F.cols(); ++i) {
			VF[F(fi, i)].push_back(fi);			
		}
	}
}

// Compute triangle-triangle adjacency with indices
template <typename DerivedF, typename DerivedTT, typename DerivedTTi>
void triangle_triangle_adjacency(
	const Eigen::PlainObjectBase<DerivedF>& F,
	Eigen::PlainObjectBase<DerivedTT>& TT,
	Eigen::PlainObjectBase<DerivedTTi>& TTi)
{
	std::vector<std::vector<int> > TTT;
//	triangle_triangle_adjacency_preprocess(F, TTT);
	for (int f = 0; f<F.rows(); ++f)
		for (int i = 0; i<F.cols(); ++i)
		{
			// v1 v2 f ei
			int v1 = F(f, i);
			int v2 = F(f, (i + 1) % F.cols());
			if (v1 > v2) std::swap(v1, v2);
			std::vector<int> r(4);
			r[0] = v1; r[1] = v2;
			r[2] = f;  r[3] = i;
			TTT.push_back(r);
		}
	std::sort(TTT.begin(), TTT.end());

//	triangle_triangle_adjacency_extractTT(F, TTT, TT);
	TT.setConstant((int)(F.rows()), F.cols(), -1);

	for (int i = 1; i<(int)TTT.size(); ++i)
	{
		std::vector<int>& r1 = TTT[i - 1];
		std::vector<int>& r2 = TTT[i];
		if ((r1[0] == r2[0]) && (r1[1] == r2[1]))
		{
			TT(r1[2], r1[3]) = r2[2];
			TT(r2[2], r2[3]) = r1[2];
		}
	}

//	triangle_triangle_adjacency_extractTTi(F, TTT, TTi);
	TTi.setConstant((int)(F.rows()), F.cols(), -1);
	for (int i = 1; i<(int)TTT.size(); ++i)
	{
		std::vector<int>& r1 = TTT[i - 1];
		std::vector<int>& r2 = TTT[i];
		if ((r1[0] == r2[0]) && (r1[1] == r2[1]))
		{
			TTi(r1[2], r1[3]) = r2[3];
			TTi(r2[2], r2[3]) = r1[3];
		}
	}
}

#ifndef IGL_STATIC_LIBRARY
//#include "fast_marching.cpp"
#endif



#endif
