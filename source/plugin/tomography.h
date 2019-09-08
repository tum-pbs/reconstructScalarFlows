/******************************************************************************
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * @author: Marie-Lena Eckert http://marielenaeckert.com/
 *
 * class for tomography solver
 *
******************************************************************************/
#ifndef TOMOGRAPHY
#define TOMOGRAPHY

#define USECG false

#include "image.h"

using namespace Eigen;

namespace Manta {
	namespace TomographyNS {

		// -------------------- helper functions grids - vecs --------------------

		KERNEL(idx) void copyGridToVec(const FlagGrid& vh, VectorX& z, const Grid<Real>& density, const int N = 0) {
			if (vh(idx) >= 0) z[N + vh(idx)] = density(idx);
		}

		KERNEL(idx) void copyVecToGridZero(Grid<Real>& density, const FlagGrid& vh, const VectorX& z, const int N = 0) {
			if (vh(idx) >= 0) density(idx) = z[N + vh(idx)];
			else density(idx) = 0;
		}

		// -------------------- tomography error calculation --------------------

		// calculate error between input image and projected density, might not be accurate as projDen only contains pixels inside the pixel visual hull
		KERNEL() void calcError(const FlagGrid* vhP, VectorX& error, const VectorX& projDen, const Grid<Real>& imgs) {
			int idx = vhP->index(i, j, k);
			if ((*vhP)(idx) == -1) error[idx] = fabs(imgs(idx)); // might be wrong if non-zero pixel in projDen outside vhP
			else error[idx] = fabs(projDen[(*vhP)(idx)] - imgs(idx));
		}

		// -------------------- copy imgs grid into vector b --------------------

		KERNEL() void setupB(const FlagGrid& vhP, VectorX& b, const Grid<Real>& imgs, const Real angleWeight) {
			if (vhP(i, j, k) != -1) {
				if (k == 0) b[vhP(i, j, k)] = imgs(i, j, k);
				else b[vhP(i, j, k)] = angleWeight*imgs(i, j, k);
			}
		}

		// -------------------- non-neg tomography --------------------

		KERNEL(idx, reduce = max) returns(Real maxVal = 0)
		Real ensureNN(const FlagGrid& vh, VectorX& z, const Grid<Real>* initialDensity) {
			if (vh(idx) >= 0) {
				Real val = z[vh(idx)];
				if (initialDensity) val += (*initialDensity)(idx);
				if (val >= 0) return;
				if (-val > maxVal) maxVal = -val;
				z[vh(idx)] = initialDensity ? -(*initialDensity)(idx) : 0;
			}
		}

		// -------------------- copy images, use front as side or back -------------------

		KERNEL() void mirrorImgToBack(Grid<Real>& imgs, const int sX, const int k0, const int k1) {
			if (k == k0 && i < sX/2) {
				imgs(i, j, k1) = imgs(sX - i - 1, j, k0);
				imgs(sX - i - 1, j, k1) = imgs(i, j, k0);
			}
		}

		KERNEL() void copyImgToSide(Grid<Real>& imgs, const int k0, const int k1) {
			if (k == k0) imgs(i, j, k1) = imgs(i, j, k0);
		}

		KERNEL() void copyImg(Grid<Real>& imgs, SimpleImage& imgPPM, const int k0) {
			if (k == k0) imgs(i, j, k) = imgPPM.get(i, j).x;
		}

	}// end namespace

	class Tomography {

	private:
		
#if USECG
		ConjugateGradient<SparseMatrix<Real>> m_lscg;
#else
		LeastSquaresConjugateGradient<SparseMatrix<Real>> m_lscg;
#endif
		
		SparseMatrix<Real, RowMajor> m_P_single; // unused if not USECG
		SparseMatrix<Real, RowMajor> m_P;
		SparseMatrix<Real, ColMajor> m_reg;
		VectorX m_b;
		VectorX m_rhsUpdate;
		VectorX m_x;
		VectorX m_z;
		VectorX m_zPrev;
		VectorX m_y;
		FlagGrid* m_vhP;
		FlagGrid* m_vh;
		const Image* m_i;

	public:

		int getVHEntry(int idx) const;

		Real computeError(const Grid<Real>& density, const Grid<Real>& imgs);

		~Tomography();

		Tomography(TomoParams& params, const Grid<Real>& imgs, const Image& i, const FlagGrid& flags, const Grid<Real>* densityMask = nullptr);

		// init before each new solve
		inline void init(const FlagGrid& flags, const Grid<Real>& imgs, const Grid<Real>* densityMask = nullptr);

		void solve(Grid<Real>& density, const Grid<Real>& imgs, const Grid<Real>* initialDen = nullptr, const FlagGrid* flags = nullptr);

		TomoParams* m_params;

	}; // class 
}

#endif