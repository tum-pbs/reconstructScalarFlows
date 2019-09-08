/******************************************************************************
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * @author: Marie-Lena Eckert http://marielenaeckert.com/
 *
 * class for divergence-free optical flow solver
 *
******************************************************************************/
#ifndef DIV_FREE_OF
#define DIV_FREE_OF

#include "pdSolve.h"
#include "tomography.h"

using namespace Eigen;

namespace Manta {
	namespace DivFreeOFNS {

		// -------------------- helper functions grids - vecs --------------------

		KERNEL(idx) void copyGridToVec(const MACGrid& velocity, const Grid<Vec3i>& vh, VectorX& z, const OFParams& ofParams) {
			for (int c = 0; c < ofParams.getDim(); c++) {
				if (vh(idx)[c] >= 0) z[vh(idx)[c]] = velocity(idx)[c];
			}
		}

		KERNEL(idx) void copyVecToGrid(MACGrid& velocity, const Grid<Vec3i>& vh, const VectorX& z, const OFParams& ofParams) {
			for (int c = 0; c < ofParams.getDim(); c++) {
				if (vh(idx)[c] >= 0) velocity(idx)[c] = z[vh(idx)[c]];
				else velocity(idx)[c] = 0;
			}
		}

		KERNEL() void setInitialDensityClearDenU(Grid<Real>& initialDensity, const Grid<Real>& denPredict, Grid<Real>& denCorrection, const FlagGrid& vh, const int NVel, const Tomography& tomo, const VectorX& z_denUpdate) {
			int idx = denPredict.index(i, j, k);
			// initialDensity is sum of denPredict and denUpdate (which values are currently in m_z of OF solve)
			if (tomo.getVHEntry(idx) >= 0) initialDensity(idx) = denPredict(idx) + z_denUpdate[NVel + vh(idx)];
			// make sure initialDensity already contains hard constraints on denUpdate (sum of predict and update is zero outside of tomo visual hull)
			else if (vh(idx) >= 0) initialDensity(idx) = 0;
			// outside any visual hull: denUpdate is zero, so initialDensity is only denPredict, should be very close to zero
			else initialDensity(idx) = 0;// denPredict(idx);

			// denCorrection is initialized to zero in order to save the denCorrection from tomography in it
			denCorrection(idx) = 0;
		}

		KERNEL() void updateDenUpdate(Grid<Real>& denCorrectionDenUpdate, VectorX& z_denUpdate, const Grid<Real>& denPredict, const FlagGrid& vh, const int NVel, const Tomography& tomo) {
			int idx = denCorrectionDenUpdate.index(i, j, k);
			if (tomo.getVHEntry(idx) >= 0) {
				// add correction onto denUpdate
				z_denUpdate[NVel + vh(idx)] += denCorrectionDenUpdate(idx); // here, denU already has the values from m_z from tomography
			}
			else if (vh(idx) >= 0) {
				// make sure sum of denUpdate and denPredict is zero outside of tomo visual hull
				z_denUpdate[NVel + vh(idx)] = -denPredict(idx);
			}

			if (vh(idx) >= 0) denCorrectionDenUpdate(idx) = z_denUpdate[NVel + vh(idx)];
			else denCorrectionDenUpdate(idx) = 0;
		}

		// -------------------- OF error calculation --------------------

		// calculate error from brightness constancy assumption / advection/transport equation
		KERNEL(bnd = 1) void calcError(const MACGrid& velocity, VectorX& error, const Grid<Real>& den, const OFParams& ofParams, const Grid<Real>* den1 = nullptr) {
			// iterate over c-component of velocity of current cell
			int idx = velocity.index(i, j, k);
			Vec3 grad, vel;
			for (int c = 0; c < ofParams.getDim(); c++) {
				grad[c] = 0.5*(den(idx + ofParams.getStride(c)) - den(idx - ofParams.getStride(c)));
				vel[c] = 0.5*(velocity(idx)[c] + velocity(idx + ofParams.getStride(c))[c]);
			}
			if (den1) error[idx] = fabs((*den1)(idx) - den(idx) + dot(grad, vel));
			else error[idx] = fabs(den(idx) + dot(grad, vel));
		}

		inline Real computeError(Real& denError, const FlagGrid& flags, const MACGrid& velUpdate, const Grid<Real>& den, const OFParams& ofParams, const Grid<Real>* den1) {
			// error from equation
			VectorX error; error.setZero(ofParams.getN());
			DivFreeOFNS::calcError(velUpdate, error, den, ofParams, den1);
			Real retError = getMaxError(den, error);
			// error in density through advection
			if (den1) {
				Grid<Real> denadv(den);
				advectSemiLagrange(&flags, &velUpdate, &denadv, 2, MACCORMACK_STRENGTH);
				denadv.sub(*den1);
				denError = denadv.getMaxAbs();
			}
			return retError;
		}

		// -------------------- setup system OF --------------------

		// adaptive penalty for kinetic energy in z direction, only applicable for single view
		inline Real getKinetic(const int c, const OFParams& ofParams, const int j, const int Y) {
			if (c < 2) return ofParams.getKinetic();
			if (!ofParams.useAdaptiveKinZ()) return ofParams.getKineticZ();
			Real half = Y / 2. - 1.;
			Real dist = fabs(half - j) / (Y / 2.);
			Real factor = 10.*(dist*dist) + 1.;
			return ofParams.getKineticZ()*factor;
		}

		// setup matrix for optical flow
		KERNEL(bnd = 1, single) void setupA(const Grid<Vec3i>& vhVel, std::vector<Triplet<Real>>& triplets, const FlagGrid& flags, const Grid<Real>& den, const Grid<Real>* den1T, const OFParams& ofParams, const FlagGrid* vhDen, const Tomography* tomo) {
			int idx = den.index(i, j, k);
			if (!flags.isFluid(idx) || flags.isSrc(i,j,k)) return;

			Real weightEntry = 1; // isOutSideVH ? 100 : 1;

			// set sigma, regularizers and identity for density in bivariate case
			int idxDen = -1; 
			if (ofParams.isBivariate()) {
				TomoParams* tomoParams = tomo->m_params;
				idxDen = (*vhDen)(idx);
				if (idxDen >= 0) {
					// set regularizers
					// iterate over neighbours in all directions for smoothness regularizer
					int mult = 0;
					for (int c = 0; c < ofParams.getDim(); c++) {
						Vec3i ijkN = Vec3i(i, j, k); ijkN[c]--;
						if ((*vhDen)(ijkN) >= 0) {
							triplets.push_back(Triplet<Real>(ofParams.getNumOfVelCmp() + idxDen, ofParams.getNumOfVelCmp() + (*vhDen)(ijkN), -(*tomoParams).getSmooth()));
							mult++;
						}
						ijkN = Vec3i(i, j, k); ijkN[c]++;
						if ((*vhDen)(ijkN) >= 0) {
							triplets.push_back(Triplet<Real>(ofParams.getNumOfVelCmp() + idxDen, ofParams.getNumOfVelCmp() + (*vhDen)(ijkN), -(*tomoParams).getSmooth()));
							mult++;
						}
					}
					triplets.push_back(Triplet<Real>(ofParams.getNumOfVelCmp() + idxDen, ofParams.getNumOfVelCmp() + idxDen, weightEntry + tomoParams->getSigma() + mult * tomoParams->getSmooth() + tomoParams->getKinetic()));
				}
			}

			if (vhVel(idx) == Vec3i(-1)) return; // cell is obstacle, can't have any density gradient, has only zero velocity components, is excluded from linear system of equations

			// set gradient
			Vec3 g;
			for (int c = 0; c < ofParams.getDim(); c++) {
				if (ofParams.shouldUseDenTarget() && den1T) g[c] = 0.25*(den(idx + ofParams.getStride(c)) - den(idx - ofParams.getStride(c)) + (*den1T)(idx + ofParams.getStride(c)) - (*den1T)(idx - ofParams.getStride(c)));
				else g[c] = 0.5*(den(idx + ofParams.getStride(c)) - den(idx - ofParams.getStride(c)));
			}

			// set grad triplets
			// iterate over c-component of velocity of current cell
			for (int c = 0; c < ofParams.getDim(); c++) {
				// index of current velocity component
				const int idxVel = vhVel(idx)[c];
				const bool isValidCmp = idxVel >= 0;
				
				// set grad triplets
				if (std::fabs(g[c]) > 1e-4) {
					// index of upper neighbouring same vel comp
					Vec3i ijkN; ijkN = Vec3i(i, j, k); ijkN[c]++;
					const int idxVelN = vhVel(ijkN)[c];
					const bool isValidCmpN = idxVelN >= 0;

					// bivariate case
					if (ofParams.isBivariate() && idxDen >= 0) {
						Real x = weightEntry*g[c] / 2.;
						if (isValidCmp) {
							triplets.push_back(Triplet<Real>(idxVel, ofParams.getNumOfVelCmp() + idxDen, x));
							triplets.push_back(Triplet<Real>(ofParams.getNumOfVelCmp() + idxDen, idxVel, x));
						}
						// upper neighbour of same component
						if (isValidCmpN) {
							triplets.push_back(Triplet<Real>(idxVelN, ofParams.getNumOfVelCmp() + idxDen, x));
							triplets.push_back(Triplet<Real>(ofParams.getNumOfVelCmp() + idxDen, idxVelN, x));
						}
					}
					
					// set grad triplets
					Real xx = weightEntry*g[c] * g[c] / 4.;
					if (std::fabs(xx) > 1e-4) {
						if (isValidCmp) triplets.push_back(Triplet<Real>(idxVel, idxVel, xx));
						// upper neighbour of same component
						if (isValidCmpN) {
							if (isValidCmp) {
								triplets.push_back(Triplet<Real>(idxVelN, idxVel, xx));
								triplets.push_back(Triplet<Real>(idxVel, idxVelN, xx));
							}
							triplets.push_back(Triplet<Real>(idxVelN, idxVelN, xx));
						}
					}

					for (int cN = 0; cN < ofParams.getDim(); cN++) {
						Real xy = weightEntry*g[c] * g[cN] / 4.;
						if (cN != c && std::fabs(xy) > 1e-4) {
							// other components of same cell
							int idxC2 = vhVel(idx)[cN];
							if (idxC2 >= 0) {
								if (isValidCmp) triplets.push_back(Triplet<Real>(idxVel, idxC2, xy));
								if (isValidCmpN) triplets.push_back(Triplet<Real>(idxVelN, idxC2, xy));
							}
							Vec3i ijkC2N = Vec3i(i, j, k); ijkC2N[cN]++;
							const int idxC2N = vhVel(ijkC2N)[cN];
							if (idxC2N >= 0) {
								if (isValidCmp) triplets.push_back(Triplet<Real>(idxVel, idxC2N, xy));
								if (isValidCmpN) triplets.push_back(Triplet<Real>(idxVelN, idxC2N, xy));
							}
						}
					}
				}

				if (!isValidCmp) continue;

				// set regularizers
				// iterate over neighbours in all directions for smoothness regularizer
				int mult = 0;
				for (int cN = 0; cN < ofParams.getDim(); cN++) {
					Vec3i ijkC2N = Vec3i(i, j, k); ijkC2N[cN]--;
					if (vhVel(ijkC2N)[c] >= 0) {
						triplets.push_back(Triplet<Real>(idxVel, vhVel(ijkC2N)[c], -ofParams.getSmooth()));
						mult++;
					}
					ijkC2N = Vec3i(i, j, k); ijkC2N[cN]++;
					if (vhVel(ijkC2N)[c] >= 0) {
						triplets.push_back(Triplet<Real>(idxVel, vhVel(ijkC2N)[c], -ofParams.getSmooth()));
						mult++;
					}
				}
				triplets.push_back(Triplet<Real>(idxVel, idxVel, ofParams.getSigma() + mult * ofParams.getSmooth() + getKinetic(c, ofParams, j, vhVel.getSizeY())));

				// add sum regularizer for z-vel component (for single view, bivariate case); only add for visible den regions
				if (ofParams.getSumZero() > 0 && ofParams.isBivariate() && (*vhDen)(idx) >= 0 && c == 2) {
					for (int k2 = 0; k2 < ofParams.getS(2); k2++) {
						if (vhVel(i, j, k2)[c] >= 0 && (*vhDen)(i, j, k2) >= 0)
							triplets.push_back(Triplet<Real>(idxVel, vhVel(i, j, k2)[c], ofParams.getSumZero()*ofParams.getSumZero()));
					}
				}
			}
		}

		// setup b for optical flow
		KERNEL(bnd = 1) void setupB(const Grid<Vec3i>& vhVel, VectorX& b, const MACGrid* velCoarse, const FlagGrid& flags, const Grid<Real>& den, const OFParams& ofParams, const Grid<Real>& den1) {
			int idx = den.index(i, j, k);
			if (!flags.isFluid(idx)) return;
			if (vhVel(idx) == Vec3i(-1)) return; // cell is obstacle, can't have any density gradient, has only zero velocity components, is excluded from linear system of equations

			// set gradient
			Vec3 g;
			for (int c = 0; c < ofParams.getDim(); c++) g[c] = 0.25*(den(idx + ofParams.getStride(c)) - den(idx - ofParams.getStride(c)) + den1(idx + ofParams.getStride(c)) - den1(idx - ofParams.getStride(c)));

			// set grad triplets, b
			// iterate over c-component of velocity of current cell
			for (int c = 0; c < ofParams.getDim(); c++) {
				// index of current velocity component
				const int idxVel = vhVel(idx)[c];
				const bool isValidCmp = idxVel >= 0;
				// index of upper neighbouring same vel comp
				Vec3i ijkN; ijkN = Vec3i(i, j, k); ijkN[c]--;

				// set grad triplets and b
				if (isValidCmp) {
					if (std::fabs(g[c]) > 1e-4) b[idxVel] += 0.5 * g[c] * (den1(idx) - den(idx));

					int idxN = den.index(ijkN);
					Real gN = 0.25*(den(idxN + ofParams.getStride(c)) - den(idxN - ofParams.getStride(c)) + den1(idxN + ofParams.getStride(c)) - den1(idxN - ofParams.getStride(c)));
					if (!flags.isObstacle(ijkN) && vhVel(ijkN) != Vec3i(-1) && std::fabs(gN) > 1e-4) 
						b[idxVel] += 0.5 * gN * (den1(idxN) - den(idxN));
				}

				// extend b for multi-scale OF
				if (velCoarse) {
					int mult = 0;
					for (int cN = 0; cN < ofParams.getDim(); cN++) {
						Vec3i ijkC2N = Vec3i(i, j, k); ijkC2N[cN]--;
						if (vhVel(ijkC2N)[c] >= 0) {
							b[idxVel] -= ofParams.getSmooth() * (*velCoarse)(ijkC2N)[c];
							mult++;
						}
						ijkC2N = Vec3i(i, j, k); ijkC2N[cN]++;
						if (vhVel(ijkC2N)[c] >= 0) {
							b[idxVel] -= ofParams.getSmooth() * (*velCoarse)(ijkC2N)[c];
							mult++;
						}
					}
					b[idxVel] += mult * ofParams.getSmooth() * (*velCoarse)(idx)[c] + ofParams.getKinetic()*(*velCoarse)(idx)[c];
				}
			}
		}

		// -------------------- hull --------------------

		// only unknown (non-zero) velocity components
		// create an index for all velocity components that can be non-zero
		KERNEL(single, reduce = +) returns(int numOfVelCmp = 0)
		int setupVelVH(Grid<Vec3i>& vh, const FlagGrid& flags, const int dim) {
			if (!flags.isFluid(i, j, k) || flags.isSrc(i, j, k) || !vh.isInBounds(Vec3i(i,j,k),2)) {
				vh(i, j, k) = Vec3i(-1);
				return;
			}
			for (int c = 0; c < dim; c++) {
				Vec3i ijkN = Vec3i(i, j, k); ijkN[c]--;
				if (!flags.isFluid(ijkN) || flags.isObstacle(ijkN)) vh(i, j, k)[c] = -1;
				else vh(i, j, k)[c] = numOfVelCmp++;
			}
		}

		// all density voxels that can take a value different from zero in the bivariate OF solve
		KERNEL(single, reduce = +) returns(int numOfDenCmp = 0)
		int setupDenVH(FlagGrid& vh, const FlagGrid& flags, const Grid<Real>& denPredict, const Tomography* tomo) { 
			if (!flags.isFluid(i, j, k) || flags.isSrc(i, j, k) || (tomo->getVHEntry(flags.index(i, j, k)) < 0 && fabs(denPredict(i, j, k)) < 1e-6) || !vh.isInBounds(Vec3i(i,j,k),2)) vh(i, j, k) = -1;
			else vh(i, j, k) = numOfDenCmp++;
		}

	}// end namespace

	class DivFreeOF {

	private:

		// optical flow on staggered grid!
		void setupSystem(const MACGrid* velCoarse, const FlagGrid& flags, const Grid<Real>& den, const Grid<Real>* den1T, const Tomography* tomo);

		inline void verifyArguments(const Grid<Real>* den1, const Tomography* tomoParams);

		int m_size;
		int m_numOfDenComp;
		ConjugateGradient<SparseMatrix<Real>, Lower | Upper> m_cg;
		SparseMatrix<Real, RowMajor> m_A;
		VectorX m_b;
		VectorX m_x;
		VectorX m_z;
		VectorX m_zPrev;
		VectorX m_y; 
		OFParams* m_params;
		Grid<Vec3i>* m_vh;
		FlagGrid* m_vhDen;
		
	public:
		Grid<Real>* m_helper;

		~DivFreeOF();

		DivFreeOF(FluidSolver* s, OFParams& ofParams, const FlagGrid& flags, const MACGrid* velCoarse, const Grid<Real>& den, Grid<Real>& pressure, const Grid<Real>* den1T, const Tomography* tomo = nullptr);

		void init(FluidSolver* s, const FlagGrid& flags, const MACGrid* velCoarse, const Grid<Real>& den, Grid<Real>& pressure, const Grid<Real>* den1T, const Tomography* tomo = nullptr);

		inline void makeDivFree(MACGrid& velUpdate, const MACGrid* velPredict, FlagGrid& flags, const FlagGrid& flagsPressure, const Grid<Real>* denPredict, const Real velInflowValue, const bool setInflow = true, const bool bivariate = true);

		int solve(MACGrid& velUpdate, const MACGrid* velPredict, FlagGrid& flags, const Grid<Real>* denPredict = nullptr, const Grid<Real>* imgs = nullptr, Tomography* tomo = nullptr, const Grid<Real>* denTarget = nullptr);

	}; // class 
}

#endif