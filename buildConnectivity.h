#ifndef BUILD_CONNECTIVITY_H_
#define BUILD_CONNECTIVITY_H_
#include <Eigen/Dense>
#include <vector>
using namespace std;

//=================build adjacent info=======================
template <typename Index, typename IndexVector>
void vertex_vertex_adjacency(const Eigen::PlainObjectBase<Index>  & F,
	std::vector<std::vector<IndexVector> >& VV) {
	VV.clear();
	VV.resize(F.maxCoeff() + 1);

	// Loop over faces
	for (int i = 0; i<F.rows(); i++) {
		// Loop over this face
		for (int j = 0; j<F.cols(); j++) {
			// Get indices of edge: s --> d
			int s = F(i, j);
			int d = F(i, (j + 1) % F.cols());
			VV.at(s).push_back(d);
			VV.at(d).push_back(s);
		}
	}
	// Remove duplicates
	for (int i = 0; i<(int)VV.size(); ++i) {
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
	for (Index fi = 0; fi<F.rows(); ++fi) {
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

#if 1
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
#endif

}

// combine the above  functions:vertex_triangle_adjacency and triangle_triangle_adjacency
template <typename DerivedF, typename VVType,typename VFType, typename DerivedTT, typename DerivedTTi>
bool buildConnectivity(const Eigen::PlainObjectBase<DerivedF>& F,
	std::vector<std::vector<VVType> >& VV,
	std::vector<std::vector<VFType> >& VT,
	Eigen::PlainObjectBase<DerivedTT>& TT,
	Eigen::PlainObjectBase<DerivedTTi>& TTi)
{
	if (F.rows() == 0) return false;//if (V.rows() == 0 || F.rows() == 0) return false;
									//V_border = igl::is_border_vertex(V, F);
	vertex_vertex_adjacency(F, VV);
	vertex_triangle_adjacency(F, VT);
	std::cout << "has build the vertex_triangle connectivity!\n";
	triangle_triangle_adjacency(F, TT, TTi);
	std::cout << "has build the triangle_triangle connectivity!\n";
	return true;
}

#endif