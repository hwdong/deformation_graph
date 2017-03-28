#include "buildConnectivity.h"
#include <iostream>

#include <Eigen/Geometry>
#include <cassert>
#include "fast_marching.h"

using namespace std;

#define FM_INLINE

#define	FM_ABS(a)       ((a) > 0 ? (a) : -(a))			//!<	Returns the absolute value a
#define FM_EPSILON 1e-9
#define FM_MIN(a,b)  std::min<FM_Float>(a,b)
#define FM_MAX(a,b)  std::max<FM_Float>(a,b)

#ifdef FM_DEBUG
#define FM_ASSERT(expr) assert(expr)
#else
#define FM_ASSERT(expr)   // if(!(expr)) cerr << "Error in file " << __FILE__ << " line " << __LINE__ << "." << endl 
#endif // FM_DEBUG


bool FastMarchingData::BuildConnectivity(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F)
{
	if (V.rows() == 0 || F.rows() == 0) return false;

	buildConnectivity(F, VV,VF, TT, TTi);
	
	return true;
}

/*------------------------------------------------------------------------------*/
// Name : ComputeUpdate_SethianMethod
/**
*  \param  d1 [FM_Float] Distance value at 1st vertex.
*  \param  d2 [FM_Float] Distance value at 2nd vertex.
*  \param  a [FM_Float] Length of the 1st edge.
*  \param  b [FM_Float] Length of the 2nd edge.
*  \param  dot [FM_Float] Value of the dot product between the 2 edges.
*  \return [FM_Float] The update value.
*  \author Gabriel Peyr?
*  \date   5-26-2003
*
*  Compute the update value using Sethian's method.
*/
/*------------------------------------------------------------------------------*/
FM_INLINE
FM_Float ComputeUpdate_SethianMethod(FM_Float d1, FM_Float d2, FM_Float a, FM_Float b, FM_Float dot, FM_Float F)
{
	FM_Float t = FM_INFINITE;

	FM_Float rCosAngle = dot;
	FM_Float rSinAngle = sqrt(1 - dot*dot);

	/* Sethian method */
	FM_Float u = d2 - d1;		// T(B)-T(A)
								//	FM_ASSERT( u>=0 );
	FM_Float f2 = a*a + b*b - 2 * a*b*rCosAngle;
	FM_Float f1 = b*u*(a*rCosAngle - b);
	FM_Float f0 = b*b*(u*u - F*F*a*a*rSinAngle*rSinAngle);

	/* discriminant of the quartic equation */
	FM_Float delta = f1*f1 - f0*f2;

	if (delta >= 0)
	{
		if (FM_ABS(f2) > FM_EPSILON)
		{
			/* there is a solution */
			t = (-f1 - sqrt(delta)) / f2;
			/* test if we must must choose the other solution */
			if (t < u ||
				b*(t - u) / t < a*rCosAngle ||
				a / rCosAngle < b*(t - u) / t)
			{
				t = (-f1 + sqrt(delta)) / f2;
			}
		}
		else
		{
			/* this is a 1st degree polynom */
			if (f1 != 0)
				t = -f0 / f1;
			else
				t = -FM_INFINITE;
		}
	}
	else
		t = -FM_INFINITE;

	/* choose the update from the 2 vertex only if upwind criterion is met */
	if (u < t &&
		a*rCosAngle < b*(t - u) / t &&
		b*(t - u) / t < a / rCosAngle)
	{
		return t + d1;
	}
	else
	{
		return FM_MIN(b*F + d1, a*F + d2);
	}
}

/*------------------------------------------------------------------------------*/
// Name : ComputeUpdate_MatrixMethod
/**
*  \param  d1 [FM_Float] Distance value at 1st vertex.
*  \param  d2 [FM_Float] Distance value at 2nd vertex.
*  \param  a [FM_Float] Length of the 1st edge.
*  \param  b [FM_Float] Length of the 2nd edge.
*  \param  dot [FM_Float] Value of the dot product between the 2 edges.
*  \return [FM_Float] The update value.
*  \author Gabriel Peyr?
*  \date   5-26-2003
*
*  Compute the update value using a change of basis method.
*/
/*------------------------------------------------------------------------------*/
FM_INLINE
FM_Float ComputeUpdate_MatrixMethod(FM_Float d1, FM_Float d2, FM_Float a, FM_Float b, FM_Float dot, FM_Float F)
{
	FM_Float t;

	/* the directional derivative is D-t*L */
	Eigen::Vector2d D = Eigen::Vector2d(d1 / b, d2 / a);
	Eigen::Vector2d L = Eigen::Vector2d(1 / b, 1 / a);

	Eigen::Vector2d QL;	//Q*L
	Eigen::Vector2d QD;	//Q*L

	FM_Float det = 1 - dot*dot;		// 1/det(Q) where Q=(P*P^T)^-1

	QD[0] = 1 / det * (D[0] - dot*D[1]);
	QD[1] = 1 / det * (-dot*D[0] + D[1]);
	QL[0] = 1 / det * (L[0] - dot*L[1]);
	QL[1] = 1 / det * (-dot*L[0] + L[1]);

	/* compute the equation 'e2*t?+ 2*e1*t + e0 = 0' */
	FM_Float e2 = QL[0] * L[0] + QL[1] * L[1];			// <L,Q*L>
	FM_Float e1 = -(QD[0] * L[0] + QD[1] * L[1]);		// -<L,Q*D>
	FM_Float e0 = QD[0] * D[0] + QD[1] * D[1] - F*F;	// <D,Q*D> - F?

	FM_Float delta = e1*e1 - e0*e2;

	if (delta >= 0)
	{
		if (FM_ABS(e2) > FM_EPSILON)
		{
			/* there is a solution */
			t = (-e1 - sqrt(delta)) / e2;
			/* upwind criterion : Q*(D-t*l)<=0, i.e. QD<=t*QL */
			if (t<FM_MIN(d1, d2) || QD[0]>t*QL[0] || QD[1] > t*QL[1])
				t = (-e1 + sqrt(delta)) / e2;	// criterion not respected: choose bigger root.
		}
		else
		{
			if (e1 != 0)
				t = -e0 / e1;
			else
				t = -FM_INFINITE;
		}
	}
	else
		t = -FM_INFINITE;
	/* choose the update from the 2 vertex only if upwind criterion is met */
	if (t >= FM_MAX(d1, d2) && QD[0] <= t*QL[0] && QD[1] <= t*QL[1])
		return t;
	else
		return FM_MIN(b*F + d1, a*F + d2);
}


int get_oppisite_f_v(const Eigen::MatrixXi &F,
	const Eigen::MatrixXi &TT,
	const Eigen::MatrixXi &TTi,
	const int f,
	const int v//,const int v1,const int v2
	,int &opposite_f,int &opposite_v) {
	int j = 0;
	int v1, v2;
	while (j < 3) {
		if (F(f, j) == v) {
			j++; if (j == 3) j = 0;
			break;
		}
		j++;
	}
	FM_ASSERT(j < 3);
	opposite_f = TT(f, j);  // cur_f is the opposite face of the face f and vert.vid
	if (opposite_f == -1) return -1;
	j = TTi(f, j);  //ej' on opposite_f
	j = (j +2) % 3;
	opposite_v = F(opposite_f, j);	
	return opposite_f;
}

/*------------------------------------------------------------------------------*/
// Name : UnfoldTriangle
/**
*  \param  CurFace [FM_GeodesicFace&] Vertex to update.
*  \param  vert [FastMarchingVertex&] Current face.
*  \param  vert1 [FastMarchingVertex&] 1st neighbor.
*  \param  vert2 [FastMarchingVertex&] 2nd neighbor.
*  \return [FastMarchingVertex*] The vertex.
*  \author Gabriel Peyr?
*  \date   5-26-2003
*
*  Find a correct vertex to update \c v.
*/
/*------------------------------------------------------------------------------*/
FM_INLINE
FastMarchingVertex * UnfoldTriangle(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
	const Eigen::MatrixXi &TT, const  Eigen::MatrixXi &TTi,
	std::vector<struct FastMarchingVertex> &vertices,
	int f, int v, int v1, int v2,
	FM_Float& dist, FM_Float& dot1, FM_Float& dot2)
{
	const Eigen::Vector3d& p = V.row(v);// vert.GetPosition();
	const Eigen::Vector3d& p1 = V.row(v1);//vert1.GetPosition();
	const Eigen::Vector3d& p2 = V.row(v2);//vert2.GetPosition();

	Eigen::Vector3d e1 = p1 - p;
	FM_Float rNorm1 = e1.norm(); //~e1
	e1.normalize(); // e1 /= rNorm1;
	Eigen::Vector3d e2 = p2 - p;
	FM_Float rNorm2 = e2.norm(); // ~e2;
	e2.normalize(); // e2 /= rNorm2;

	FM_Float dot = e1.adjoint()*e2;// e1*e2;
	FM_ASSERT(dot < 0);

	/* the equation of the lines defining the unfolding region [e.g. line 1 : {x ; <x,eq1>=0} ]*/
	Eigen::Vector2d eq1 = Eigen::Vector2d(dot, sqrt(1 - dot*dot));
	Eigen::Vector2d eq2 = Eigen::Vector2d(1, 0);

	/* position of the 2 points on the unfolding plane */
	Eigen::Vector2d x1(rNorm1, 0);
	Eigen::Vector2d x2 = eq1*rNorm2;

	/* keep track of the starting point */
	Eigen::Vector2d xstart1 = x1;
	Eigen::Vector2d xstart2 = x2;

	FastMarchingVertex* pV1 = &(vertices[v1]);
	FastMarchingVertex* pV2 = &(vertices[v2]);

	
	int cur_f, cur_v;
	
	get_oppisite_f_v(F, TT,TTi, f,v, cur_f, cur_v);

	//FM_GeodesicFace* pCurFace = (FM_GeodesicFace*)CurFace.GetFaceNeighbor(vert);


	int nNum = 0;
	while (nNum < 50 && cur_f != -1) // NULL)
	{
		//	FastMarchingVertex* pV = (FastMarchingVertex*)pCurFace->GetVertex(*pV1, *pV2); //opposite vertex to face and edge(pV1,pV2)
		//	FM_ASSERT(pV != NULL);
		FastMarchingVertex* pV = &(vertices[cur_v]); //opposite vertex to face and vert
													 /*
													 e1 = pV2->GetPosition() - pV1->GetPosition();
													 FM_Float rNorm1 = ~e1;
													 e1 /= rNorm1;
													 e2 = pV->GetPosition() - pV1->GetPosition();
													 FM_Float rNorm2 = ~e2;
													 e2 /= rNorm2;
													 */

		Eigen::Vector3d e1 = V.row(pV2->vid) - V.row(pV1->vid);
		FM_Float rNorm1 = e1.norm(); //~e1
		e1.normalize(); // e1 /= rNorm1;
		Eigen::Vector3d e2 = V.row(pV->vid) - V.row(pV1->vid);
		FM_Float rNorm2 = e2.norm(); // ~e2;
		e2.normalize(); // e2 /= rNorm2;

						/* compute the position of the new point x on the unfolding plane (via a rotation of -alpha on (x2-x1)/rNorm1 )
						| cos(alpha) sin(alpha)|
						x = |-sin(alpha) cos(alpha)| * [x2-x1]*rNorm2/rNorm1 + x1   where cos(alpha)=dot
						*/
		Eigen::Vector2d vv = (x2 - x1)*rNorm2 / rNorm1;
		dot = e1.adjoint()*e2;  //e1*e2;
								//	Eigen::Vector2d x = vv.Rotate2D();////vv.Rotate(-acos(dot)) + x1;

		Eigen::Rotation2D<double> rot2(-acos(dot));
		Eigen::Vector2d x = rot2*vv + x1;  //dhw to check


										   /* compute the intersection points.
										   We look for x=x1+lambda*(x-x1) or x=x2+lambda*(x-x2) with <x,eqi>=0, so */
		FM_Float lambda11 = -(x1.dot(eq1)) / ((x - x1).dot(eq1));	 //-(x1*eq1) / ((x - x1)*eq1);	// left most 
		FM_Float lambda12 = -(x1.dot(eq2)) / ((x - x1).dot(eq2));  //-(x1*eq2) / ((x - x1)*eq2);	// right most
		FM_Float lambda21 = -(x2.dot(eq1)) / ((x - x2).dot(eq1)); //-(x2*eq1) / ((x - x2)*eq1);	// left most 
		FM_Float lambda22 = -(x2.dot(eq2)) / ((x - x2).dot(eq2));   //-(x2*eq2) / ((x - x2)*eq2);	// right most
		bool bIntersect11 = (lambda11 >= 0) && (lambda11 <= 1);
		bool bIntersect12 = (lambda12 >= 0) && (lambda12 <= 1);
		bool bIntersect21 = (lambda21 >= 0) && (lambda21 <= 1);
		bool bIntersect22 = (lambda22 >= 0) && (lambda22 <= 1);
		if (bIntersect11 && bIntersect12)
		{
			//			FM_ASSERT( !bIntersect21 && !bIntersect22 );
			/* we should unfold on edge [x x1] */
			//	pCurFace = (FM_GeodesicFace*)pCurFace->GetFaceNeighbor(*pV2);
			f = cur_f;
			get_oppisite_f_v(F, TT,TTi, f, pV2->vid, cur_f, cur_v);

			pV2 = pV;
			x2 = x;
		}
		else if (bIntersect21 && bIntersect22)
		{
			//			FM_ASSERT( !bIntersect11 && !bIntersect12 );
			/* we should unfold on edge [x x2] */
			//	pCurFace = (FM_GeodesicFace*)pCurFace->GetFaceNeighbor(*pV1);
			f = cur_f;
			get_oppisite_f_v(F, TT,TTi, f, pV1->vid, cur_f, cur_v);

			pV1 = pV;
			x1 = x;
		}
		else
		{
			FM_ASSERT(bIntersect11 && !bIntersect12 &&
				!bIntersect21 && bIntersect22);
			/* that's it, we have found the point */
			dist = x.norm(); // ~x;
			dot1 = x.dot(xstart1) / (dist * xstart1.norm());//  ~xstart1);
			dot2 = x.dot(xstart2) / (dist * xstart2.norm());// ~xstart2);
			return pV;
		}
		nNum++;
	}

	return NULL;
}

double ComputeVertexDistance(
	const Eigen::MatrixXd &V, const Eigen::MatrixXi &Faces,
	const Eigen::MatrixXi &TT, const  Eigen::MatrixXi &TTi,
	std::vector<struct FastMarchingVertex> &vertices,
	int &cur_f,
	FastMarchingVertex& CurrentVertex,
	FastMarchingVertex& Vert1, FastMarchingVertex& Vert2, const int CurrentFront,
	double F//weight of CurrentVertex
	, bool bUseUnfolding_ = true
	)
{

	if (Vert1.state != VertexState::kFar || Vert2.state != VertexState::kFar)
	{
		Eigen::Vector3d Edge1 = V.row(Vert1.vid) - V.row(CurrentVertex.vid);
		double b = Edge1.norm();		
		Edge1.normalize();
		//	Eigen::Vector3d Edge1 = Vert1.GetPosition() - CurrentVertex.GetPosition();
		//	FM_Float b = Edge1.Norm();
		//	Edge1 /= b;
		Eigen::Vector3d Edge2 = V.row(Vert2.vid) - V.row(CurrentVertex.vid);
		double a = Edge2.norm();
		Edge2.normalize(); 
		//	Eigen::Vector3d Edge2 = Vert2.GetPosition() - CurrentVertex.GetPosition();
		//	FM_Float a = Edge2.Norm();
		//	Edge2 /= a;

		double d1 = Vert1.distance;
		double d2 = Vert2.distance;

		/*	Set it if you want only to take in acount dead vertex
		during the update step. */
		// #define USING_ONLY_DEAD


#ifndef USING_ONLY_DEAD 
		bool bVert1Usable = Vert1.state != VertexState::kFar && Vert1.front == CurrentFront;
		bool bVert2Usable = Vert2.state != VertexState::kFar && Vert2.front == CurrentFront;
		if (!bVert1Usable && bVert2Usable)
		{
			/* only one point is a contributor */
			return d2 + a * F;
		}
		if (bVert1Usable && !bVert2Usable)
		{
			/* only one point is a contributor */
			return d1 + b * F;
		}
		if (bVert1Usable && bVert2Usable)
		{
#else
		bool bVert1Usable = Vert1.GetState() == FastMarchingVertex::kDead && Vert1.GetFront() == &CurrentFront;
		bool bVert2Usable = Vert2.GetState() == FastMarchingVertex::kDead && Vert2.GetFront() == &CurrentFront;
		if (!bVert1Usable && bVert2Usable)
		{
			/* only one point is a contributor */
			return d2 + a * F;
		}
		if (bVert1Usable && !bVert2Usable)
		{
			/* only one point is a contributor */
			return d1 + b * F;
		}
		if (bVert1Usable && bVert2Usable)
		{
#endif	// USING_ONLY_DEAD
			double dot = Edge1.dot(Edge2);

			/*	you can choose wether to use Sethian or my own derivation of the equation.
			Basicaly, it gives the same answer up to normalization constants */
#define USE_SETHIAN

			/* first special case for obtuse angles */
			if (dot < 0 && bUseUnfolding_)
			{
				double c, dot1, dot2;
				FastMarchingVertex* pVert = UnfoldTriangle(V, Faces,TT, TTi, vertices, cur_f, CurrentVertex.vid, Vert1.vid, Vert2.vid, c, dot1, dot2);
				//	FastMarchingVertex* pVert = UnfoldTriangle(CurrentFace, CurrentVertex, Vert1, Vert2, c, dot1, dot2);
				if (pVert != NULL && pVert->state != VertexState::kFar)
				{
					double d3 = pVert->distance;
					double t;		// newly computed value
									/* use the unfolded value */
#ifdef USE_SETHIAN
					t = ComputeUpdate_SethianMethod(d1, d3, c, b, dot1, F);
					t = std::min<double>(t, ComputeUpdate_SethianMethod(d3, d2, a, c, dot2, F));
#else
					t = ComputeUpdate_MatrixMethod(d1, d3, c, b, dot1, F);
					t = FM_MIN(t, ComputeUpdate_MatrixMethod(d3, d2, a, c, dot2, F));
#endif
					return t;
				}
			}

#ifdef USE_SETHIAN
			return ComputeUpdate_SethianMethod(d1, d2, a, b, dot, F);
#else
			return ComputeUpdate_MatrixMethod(d1, d2, a, b, dot, F);
#endif

		}
	}

	return FM_INFINITE;
}
	bool FastMarchingData::frontOverlap(FastMarchingVertex* pNewVert, 
		int new_front, FM_Float new_distance)
	{		
		if (pNewVert->front == new_front) return false;		
		typedef  FastMarchingVertex::FrontOverlapInfo FrontOverlapInfo;		
		if (!pNewVert->overlapInfo) {
			pNewVert->overlapInfo.reset(new FrontOverlapInfo());
			pNewVert->overlapInfo->fronts.push_back(
				FastMarchingVertex::Front_Distance(pNewVert->front, pNewVert->distance));
			pNewVert->overlapInfo->fronts.push_back(
				FastMarchingVertex::Front_Distance(new_front, new_distance));
		}
		else {
			int i = 0;
			for (; i<pNewVert->overlapInfo->fronts.size(); i++) {
				if (pNewVert->overlapInfo->fronts[i].i == new_front&&pNewVert->overlapInfo->fronts[i].value>new_distance)
					pNewVert->overlapInfo->fronts[i].value = new_distance;				
			}
			if (i == pNewVert->overlapInfo->fronts.size())
				pNewVert->overlapInfo->fronts.push_back(
					FastMarchingVertex::Front_Distance(new_front, new_distance));
		}
		return true;
	}

	bool FastMarchingData::toStop(FastMarchingVertex* pVert) {
		if (pVert->distance > option.distmax) return true;
		for (int k = 0; k < option._end_points.size(); ++k)
			if (option._end_points[k] == pVert->vid)
				return true;
		return false;
	}

bool FastMarchingData::PrepareFastMarching(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F) 
{
	int nV = V.rows();

	if (TT.rows() != nV)
		BuildConnectivity(V, F);

	if (vertices.size() > 0) vertices.clear();
	vertices.resize(nV);
	for (int i = 0; i < nV; i++) {
		vertices[i].vid = i;	
		vertices[i].reset(); //	vertices[i].front = -1;
	}

	seed_points.clear();
	return true;
}

bool FastMarchingData::PerformFastMarching(
	const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
	const vector<int> &start_points	)
{
	nbr_iter = 0;
	for (int i = 0; i < V.rows(); i++) {
		vertices[i].state = VertexState::kFar;
	}

	//set start points	
	VertexPtrVector active_Vertices;
	for (auto i : start_points) {
		vertices[i].state = VertexState::kAlive;
		vertices[i].front = i;
		vertices[i].distance = 0;
		active_Vertices.push_back(&vertices[i]);

		seed_points.push_back(i);
	}

	std::make_heap(active_Vertices.begin(), active_Vertices.end(),
		[](FastMarchingVertex* vetex1, FastMarchingVertex* vetex2)->bool 
	           {	return vetex1->distance > vetex2->distance; });

	int iter = 0;
	while (!active_Vertices.empty() && iter < option.iter_max) {
		FastMarchingVertex* pCurVert = active_Vertices.front();
		assert(pCurVert != 0);
		std::pop_heap(active_Vertices.begin(), active_Vertices.end(),
			[](FastMarchingVertex* vetex1, FastMarchingVertex* vetex2)->bool { 
			return vetex1->distance > vetex2->distance; });
		active_Vertices.pop_back();
		pCurVert->state = VertexState::kDead;

		if (pCurVert->distance >= option.distmax) 		
			return false;   //DHW_
		NewDeadVertexCallback_(pCurVert);

		// update front
		for (int i : VV[pCurVert->vid]) {
			FastMarchingVertex* pNewVert = &(vertices[i]);
			if (pCurVert->stopped&&!pNewVert->stopped&&pNewVert->state == VertexState::kFar)
				continue;
	
			/* compute it's new distance using neighborhood information */
			double new_distance = FM_INFINITE;
			int j, k;
			for (int f : VF[pNewVert->vid]) {
				if (pNewVert->vid == F(f, 0)) { j = F(f, 1); k = F(f, 2); }
				else if (pNewVert->vid == F(f, 1)) { j = F(f, 2); k = F(f, 0); }
				else { j = F(f, 0); k = F(f, 1); }
				FastMarchingVertex *pVert1 = &(vertices[j]);
				FastMarchingVertex *pVert2 = &(vertices[k]);
				assert(pVert1 != NULL&&pVert2 != NULL);
				if (pVert1->distance > pVert2->distance) {
					FastMarchingVertex* pTempVert = pVert1;
					pVert1 = pVert2;
					pVert2 = pTempVert;
				}
				double w = 1.;
				if (option.weight != 0) w = option.weight[i];
				if(pVert1->distance<new_distance)
				new_distance = std::min<double>(new_distance,
					ComputeVertexDistance(V, F, TT, TTi, vertices, f, *pNewVert, *pVert1, *pVert2, pCurVert->front, w));
				//ComputeVertexDistance(*pFace, *pNewVert, *pVert1, *pVert2, *pCurVert->GetFront()));
			}


			switch (pNewVert->state) {
			case VertexState::kFar:
				/* ask to the callback if we should update this vertex and add it to the path */
		//		if (VertexInsersionCallback_ == NULL ||	VertexInsersionCallback_(*pNewVert, new_distance))
				if(VertexInsersionCallback_(pNewVert, new_distance))
				{
				//	pNewVert->SetDistance(rNewDistance);
					pNewVert->distance = new_distance;

					/* add the vertex to the heap */
					active_Vertices.push_back(pNewVert);
					std::push_heap(active_Vertices.begin(), active_Vertices.end(), [](FastMarchingVertex* vetex1,
						FastMarchingVertex* vetex2)->bool { return vetex1->distance > vetex2->distance; });
					/* this one can be added to the heap */
					pNewVert->state = VertexState::kAlive;
					pNewVert->front = pCurVert->front;
				}
				break;
			case VertexState::kAlive:
				// just update it's value 
				if (new_distance <= pNewVert->distance)
				{
					// possible overlap with old value 
			//		if (pCurVert->front != pNewVert->front)->GetFrontOverlapInfo().RecordOverlap(*pNewVert->GetFront(), pNewVert->GetDistance());
					frontOverlap(pNewVert, pCurVert->front, new_distance);
					pNewVert->distance = new_distance;
					pNewVert->front = pCurVert->front;
					// hum, check if we can correct this (avoid recomputing the whole heap).
					std::make_heap(active_Vertices.begin(), active_Vertices.end(), [](FastMarchingVertex* vetex1,
						FastMarchingVertex* vetex2)->bool { return vetex1->distance > vetex2->distance; });
				}
				else
				{
					// possible overlap with new value 
				//	if (pCurVert->front != pNewVert->front)pNewVert->GetFrontOverlapInfo().RecordOverlap(*pCurVert->GetFront(), rNewDistance);
					frontOverlap(pNewVert, pCurVert->front, new_distance);
				}
				break;
			case VertexState::kDead:
				/* inform the user if there is an overlap */
				//if (pCurVert->front != pNewVert->front)pNewVert->GetFrontOverlapInfo().RecordOverlap(*pCurVert->GetFront(), rNewDistance);
				frontOverlap(pNewVert, pCurVert->front, new_distance);
				break;
			default:
				FM_ASSERT(false);
			}
		//	if (new_distance >= option.distmax) 			return false;   //DHW_
		}
	}
	return false;	
}

std::pair<int, FM_Float> FastMarchingData::FarthestPointSamplingStep(
	const Eigen::MatrixXd &V,
	const Eigen::MatrixXi &F) 
{
	if (!option.bound)option.bound = new double[V.rows()];

	FM_Float mad_dist = 0; int max_ind = 0;
	for (int j = 0; j < V.rows(); j++) {
		option.bound[j] = vertices[j].distance;
		if (vertices[j].distance > mad_dist) {
			mad_dist = vertices[j].distance;
			max_ind = j;
		}
	}
	vector<int> startPoints; startPoints.push_back(max_ind);
	PerformFastMarching(V, F, startPoints);
	std::cout << "max dist:=" << mad_dist << "\n";

	return std::make_pair(max_ind, mad_dist); 
}

std::pair<int, FM_Float> FastMarchingData::FarthestPointSampling(
	const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
	const vector<int> &start_points, int num) {
	PerformFastMarching(V, F, start_points);
	
	std::pair<int, FM_Float> p;
	for (int i = 0; i < num; i++){
		p  = FarthestPointSamplingStep(V, F);
	}
	return p;
}

