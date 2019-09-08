/******************************************************************************
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * @author: Marie-Lena Eckert http://marielenaeckert.com/
 *
 * Reconstructions of density, velocity, and inflow -- main methods
 *
 ******************************************************************************/

#include "tomography.h"
#include "divFreeOF.h"
#include "particle.h"

using namespace std;

namespace Manta {

	// -------------------- externally defined functions --------------------

	extern void interpolateGrid(Grid<Real>& target, const Grid<Real>& source, Vec3 scale = Vec3(1.), Vec3 offset = Vec3(0.), Vec3i size = Vec3i(-1, -1, -1), int orderSpace = 1);
	extern void interpolateMACGrid(MACGrid& target, const MACGrid& source, Vec3 scale = Vec3(1.), Vec3 offset = Vec3(0.), Vec3i size = Vec3i(-1, -1, -1), int orderSpace = 1);
	extern void cgSolveDiffusion(const FlagGrid& flags, GridBase& grid, Real alpha = 0.25, Real cgMaxIterFac = 1.0, Real cgAccuracy = 1e-4);
	extern void setWallBcs(const FlagGrid& flags, MACGrid& vel, const MACGrid* obvel = 0, const MACGrid* fractions = 0, const Grid<Real>* phiObs = 0, int boundaryWidth = 0);
	extern int blurRealGrid(const Grid<Real>& oG, Grid<Real>& tG, float si);

	// --------------------------------------------------------------------------------------
	// -------------------- combined density and velocity reconstruction --------------------
	// --------------------------------------------------------------------------------------

	// -------------------- loading and assembling images for tomography --------------------

	PYTHON() void loadPPMToGrid(Grid<Real>& imgs, const string filenamePPM, const int k0) {
		if (k0 > imgs.getSizeZ() - 1) errMsg("Error in loadPPMToGrid: #images=" << imgs.getSizeZ() << ", k0=" << k0);
		SimpleImage imgPPM = SimpleImage();
		imgPPM.init(imgs.getSizeX(), imgs.getSizeY());
		imgPPM.initFromPpm(filenamePPM);
		TomographyNS::copyImg(imgs, imgPPM, k0);
	}

	PYTHON() void reuseFrontImg(Grid<Real>& imgs, const string filenamePPM, const Image& i, const bool asSide, const bool asBack) {
		// error checking
		if (i.getCamCount() != 1 + asSide + asBack) errMsg("Error in reuseFrontImg: #cameras=" << i.getCamCount() << ", asSide=" << asSide << ", asBack=" << asBack);

		// load front image
		loadPPMToGrid(imgs, filenamePPM, 0);
		// copy images
		if (asSide) TomographyNS::copyImgToSide(imgs, 1, 0);
		if (asBack && asSide) TomographyNS::mirrorImgToBack(imgs, imgs.getSizeX(), 2, 0);
		if (asBack && !asSide) TomographyNS::mirrorImgToBack(imgs, imgs.getSizeX(), 1, 0);
	}

	// -------------------- general helper functions --------------------

	KERNEL(bnd=1) void setVelInflowKN(MACGrid& vel, const FlagGrid& flags, const Vec3& value) {
		if (flags.isSrc(i, j, k)) vel(i, j, k) = value;
		else if (flags.isSrc(i, j - 1, k)) vel(i, j - 1, k).y = value.y;
	}
	PYTHON() void setVelInflow(MACGrid& vel, const FlagGrid& flags, const Vec3& value) {
		setVelInflowKN(vel, flags, value);
	}

	KERNEL() void deleteSmokeKn(Grid<Real>& den, const Vec3i& p0, const Vec3i& p1) {
		int off = 0.05*den.getSizeX();
		// further away from source
		if (j >= p1.y && j < p1.y + 0.02*den.getSizeY() && (i<p0.x - off || i>p1.x + off || k<p0.z - off || k>p1.z + off)) den(i, j, k) = 0;

		off = 0.06*den.getSizeX();
		// further away from source
		if (j >= p1.y && j < p1.y + 0.03*den.getSizeY() && (i<p0.x - off || i>p1.x + off || k<p0.z - off || k>p1.z + off)) den(i, j, k) = 0;

		off = 0.08*den.getSizeX();
		// further away from source
		if (j >= p1.y && j < p1.y + 0.04*den.getSizeY() && (i<p0.x - off || i>p1.x + off)) den(i, j, k) = 0;

		off = 0.15*den.getSizeX();
		// further away from source
		if (j >= p1.y && j < p1.y + 0.08*den.getSizeY() && (i<p0.x - off || i>p1.x + off)) den(i, j, k) = 0;

		// directly at source
		else if (j >= p1.y && j < p1.y + 1 && (i<p0.x || i>p1.x || k<p0.z || k>p1.z)) den(i, j, k) = 0;


		if (j < p1.y) den(i, j, k) = 0;
	}

	PYTHON() void deleteSmoke(Grid<Real>& den, const Vec3i& p0, const Vec3i& p1) {
		deleteSmokeKn(den, p0, p1);
	}

	KERNEL() void deleteInflowDenKn(Grid<Real>& den, const Vec3i& p0, const Vec3i& p1, const int plus) {
		if (j <= p1.y+plus) den(i, j, k) = 0;
	}

	PYTHON() void deleteInflowDen(Grid<Real>& den, const Vec3i& p0, const Vec3i& p1, const int plus=1) {
		deleteInflowDenKn(den, p0, p1, plus);
	}

	KERNEL() void deleteInflowVelKn(MACGrid& vel, const Vec3i& p0, const Vec3i& p1, const int plus) {
		if (j <= p1.y+plus) vel(i, j, k) = Vec3(0.);
	}

	PYTHON() void deleteInflowVel(MACGrid& vel, const Vec3i& p0, const Vec3i& p1, const int plus=1) {
		deleteInflowVelKn(vel, p0, p1, plus);
	}

	KERNEL(idx, reduce = max) returns(Real maxVal = -std::numeric_limits<Real>::max())
		Real CompMaxReal(const Grid<Real>& val) {
		if (val[idx] > maxVal)
			maxVal = val[idx];
	}

	KERNEL(reduce = max) returns(Real maxVal = -std::numeric_limits<Real>::max())
	Real deleteInflowImgKn(const Grid<Real>& imgs0, const Grid<Real>& imgs1, Grid<Real>& imgs2, const int kCurrent) {
		if (k==kCurrent && fabs(imgs0(i,j,k) - imgs1(i, j, k)) > 1e-9) {
			imgs2(i, j, k) = 0;
			if (j > maxVal) maxVal = j;
		}
	}

	KERNEL() void deleteBelow(Grid<Real>& imgs2, const Real highestJ, const int kCurrent) {
		if(j<=highestJ && k == kCurrent) imgs2(i,j,k) = 0;
	}

	PYTHON() void deleteInflowImg(const Grid<Real>& imgs0, const Grid<Real>& imgs1, Grid<Real>& imgs2) {
		for (int k = 0; k < imgs0.getSizeZ(); k++) {
			Real highestJ = deleteInflowImgKn(imgs0, imgs1, imgs2, k);
			deleteBelow(imgs2, highestJ, k);
		}
	}

	// add source to density grid
	KERNEL(idx) void knAddSrc(Grid<Real>& density, const Grid<Real>& src) {
		if (density(idx) < src(idx)) density(idx) = src(idx);
	}

	// add source to density grid
	PYTHON() void addSrc(Grid<Real>& density, const Grid<Real>& src) {
		knAddSrc(density, src);
	}

	KERNEL() void deleteSmokeInBndsAdaptFlagsKN(Grid<Real>& den, FlagGrid& flags, const int bWidth = 1) {
		if (flags.isObstacle(i, j, k) || !flags.isInBounds(Vec3i(i, j, k), bWidth + 1) || flags.isOutflow(i, j, k)) den(i, j, k) = 0;
		if (!flags.isInBounds(Vec3i(i, j, k), bWidth + 1) || flags.isOutflow(i, j, k)) {
			flags(i, j, k) = FlagGrid::TypeEmpty | FlagGrid::TypeOutflow;
		}
	}

	KERNEL() void deleteSmokeInBndsKN(Grid<Real>& den, const FlagGrid& flags, const int bWidth = 1) {
		if (flags.isObstacle(i, j, k) || !flags.isInBounds(Vec3i(i, j, k), bWidth + 1) || flags.isOutflow(i, j, k)) den(i, j, k) = 0;
	}

	KERNEL() void upsampleFlagGridKN(FlagGrid& target, const FlagGrid& source, const int factor) {
		target(i, j, k) = source(int(floor(i / factor)), int(floor(j / factor)), int(floor(k / factor)));
	}

	PYTHON() void upsampleFlagGrid(FlagGrid& target, const FlagGrid& source) {
		int factor = target.getSizeX() / source.getSizeX();
		if (target.getSizeX() < source.getSizeX()) errMsg("Error in upsampleFlagGrid: resolution of target is smaller than resolution of source.");
		upsampleFlagGridKN(target, source, factor);
	}

	KERNEL()
	void downsampleFlagGrid(FlagGrid& target, const FlagGrid& source, int orderSpace = 1) {
		Vec3 sourceFactor = calcGridSizeFactor(source.getSize(), target.getSize());
		Vec3 pos = Vec3(i, j, k) * sourceFactor + sourceFactor*0.5;
		if (!source.is3D()) pos[2] = 0; // allow 2d -> 3d

		// attention! source factor is assumed to be 2 here.. 
		Vec3i neighbors[8];
		neighbors[0] = Vec3i(2 * i, 2 * j, 2 * k);
		neighbors[1] = Vec3i(2 * i + 1, 2 * j, 2 * k);
		neighbors[2] = Vec3i(2 * i, 2 * j + 1, 2 * k);
		neighbors[3] = Vec3i(2 * i, 2 * j, 2 * k + 1);
		neighbors[4] = Vec3i(2 * i + 1, 2 * j + 1, 2 * k);
		neighbors[5] = Vec3i(2 * i + 1, 2 * j, 2 * k + 1);
		neighbors[6] = Vec3i(2 * i, 2 * j + 1, 2 * k + 1);
		neighbors[7] = Vec3i(2 * i + 1, 2 * j + 1, 2 * k + 1);

		bool isOut = false;
		bool isSrc = false;
		bool isInfl = false;
		for (int n = 0; n < 8; n++) {
			// obstacle is "strongest flag"
			if (source.isObstacle(neighbors[n])) {
				target(i, j, k) = FlagGrid::TypeObstacle;
				return;
			}
			if (source.isOutflow(neighbors[n])) isOut = true;
			if (source.isSrc(neighbors[n])) isSrc = true;
			if (source.isInflow(neighbors[n])) isInfl = true;
		}

		// ordering is important here!
		if (isOut) {
			target(i, j, k) = FlagGrid::TypeOutflow | FlagGrid::TypeEmpty;
			return;
		}
		if (isInfl) {
			target(i, j, k) = FlagGrid::TypeInflow;
			return;
		}
		if (isSrc) {
			target(i, j, k) = FlagGrid::TypeSrc;
			return;
		}

		target(i, j, k) = source.getInterpolatedHi(pos, orderSpace);
		// check if target isn't a flag or a combination of those
		if (!(target.isFluid(i, j, k) || target.isObstacle(i, j, k) || target.isEmpty(i, j, k) || target.isOutflow(i, j, k) || target.isSrc(i, j, k))) {
			errMsg("Error in downsampleFlagGrid: target grid is not a valid (used) type "<<target(i,j,k));
		}
		if (!target.isInBounds(Vec3i(i, j, k), 1) && !(target.isObstacle(i, j, k) || target.isOutflow(i, j, k))) errMsg("Error in downsampleFlagGrid: bnd is neither obstacle nor outflow");
	}

	KERNEL()
	void interpolateImgsKN(Grid<Real>& target, Grid<Real>& source, const Real scale, int orderSpace = 1) {
		if (i == 0 || j == 0 || i == target.getSizeX() - 1 || j == target.getSizeY() - 1) target(i, j, k) = 0;
		else {
			Vec3 sourceFactor = calcGridSizeFactor(source.getSize(), target.getSize());
			Vec3 pos = Vec3(i, j, k) * sourceFactor + sourceFactor*0.5;
			pos.z = k + 0.5;
			target(i, j, k) = scale * source.getInterpolatedHi(pos, orderSpace);
		}
	}

	PYTHON() void interpolateImgs(Grid<Real>& target, Grid<Real>& source, const Real scale) {
		interpolateImgsKN(target, source, scale);
	}

	// -------------------- set inflow area --------------------

	KERNEL(bnd=2) void setInflowStructureKN(FlagGrid& flags, const Vec3& p0, const Vec3& p1, const Real& value, const bool noise, Grid<Real>* den) {
		if ((i <= p0.x || i >= p1.x || k <= p0.z || k >= p1.z) && j <= p1.y) {
			flags(i, j, k) = FlagGrid::TypeObstacle;
		}
		else if (j <= p1.y && !flags.isOutflow(i, j, k)) {
			flags(i, j, k)      = FlagGrid::TypeSrc | FlagGrid::TypeFluid;
			if (den) (*den)(i, j, k) = value + (noise ? rand() / (0.5*RAND_MAX) - 1. : 0);
		}
	}

	PYTHON() void setInflowStructure(FlagGrid& flags, const ShapeDetails& srcInflow, const Real& value = 10, const bool noise = false, Grid<Real>* den = nullptr) {
		srand(0);
		if (srcInflow.shape != 2) errMsg("Error in setInflowStructure: only implemented for box.");

		Vec3 p0 = srcInflow.center - srcInflow.vec;
		Vec3 p1 = srcInflow.center + srcInflow.vec;
		setInflowStructureKN(flags, p0, p1, value, noise, den);
	}

	KERNEL(idx, reduce = +) returns(Real sum = 0)
	Real mseKN_vel(const MACGrid& velDiff, const int dim, const Grid<Real>* density) {
		if (!density || fabs((*density)(idx)) > 1e-9) {
			for (int c = 0; c < dim; c++)
				sum += pow(velDiff(idx)[c], 2);
		}
	}

	KERNEL(idx, reduce = +) returns(Real sumCells = 0)
	Real cntCells(const Grid<Real>& density) {
		if (fabs(density(idx)) > 1e-9) sumCells++;
	}

	PYTHON() Real psnr_vel(const MACGrid& velDiff, const Real mx = 10, const Grid<Real>* density = 0) {
		Real scale = density ? 1. / cntCells(*density) : (1. / (velDiff.getSizeX()*velDiff.getSizeY()*velDiff.getSizeZ()));
		return 20 * log10(mx) - 10 * log10(scale*mseKN_vel(velDiff, velDiff.is3D() ? 3 : 2, density));
	}

	KERNEL(idx, reduce = +) returns(Real sum = 0)
	Real mseKN_den(const Grid<Real>& imgDiff) {
		sum += pow(imgDiff(idx), 2);
	}

	PYTHON() Real psnr_den(const Grid<Real>& imgDiff, const Real mx = 10) {
		Real scale = (1. / (imgDiff.getSizeX()*imgDiff.getSizeY()*imgDiff.getSizeZ()));
		return 20 * log10(mx) - 10 * log10(scale*mseKN_den(imgDiff));
	}

	KERNEL(idx, reduce=+) returns(Real sum = 0)
	Real l2_denKN(const Grid<Real>& den) {
		sum += den(idx) * den(idx);
	}
	PYTHON() Real l2_den(const Grid<Real>& den) {
		return sqrt(l2_denKN(den)) / (den.getSizeX()*den.getSizeY()*den.getSizeZ());
	}

	KERNEL(idx, reduce=+) returns (Real sum = 0)
	Real l2_velKN(const MACGrid& vel) {
		sum += (vel(idx).x*vel(idx).x + vel(idx).y*vel(idx).y + vel(idx).z*vel(idx).z);
	}

	PYTHON() Real l2_vel(const MACGrid& vel) {
		return sqrt(l2_velKN(vel));// / (vel.getSizeX()*vel.getSizeY()*vel.getSizeZ());
	}

	// -------------------- prediction and alignment --------------------
	
	KERNEL(idx) void setSrcToConst(Grid<Real>& src0, const FlagGrid& flags, const Real value) {
		if (flags.isSrc(idx)) src0(idx) = value;
	}

	KERNEL(bnd = 2) void delOutsideSrc(Grid<Real>& den, const FlagGrid& vh, const FlagGrid& flags) {
		if (flags.isSrc(i, j, k) && vh(i, j, k) < 0) den(i, j, k) = 0;
	}

	KERNEL(bnd = 2) void updateSrc(Grid<Real>& den, const FlagGrid& vh, const FlagGrid& flags, VectorX& x) {
		if (flags.isSrc(i, j, k) && vh(i, j, k) >= 0) den(i, j, k) += x[vh(i, j, k)];
	}

	KERNEL(bnd = 2) void extrapolateSrc(Grid<Real>& den, const FlagGrid& vh, const FlagGrid& flags) {
		if (flags.isSrc(i, j, k) && vh(i, j, k) == -1){
			if (flags.isSrc(i, j + 1, k) && vh(i, j + 1, k) != -2) {
				Real sum = 0;
				Real num = 0;
				for (int z = -1; z <= 1; z++) {
					for (int y = -1; y <= 1; y++) {
						for (int x = -1; x <= 1; x++) {
							if (vh(i + x, j + y, k + z) >= 0) {
								sum += den(i + x, j + y, k + z);
								num++;
							}
						}
					}
				}
				if (num>0) den(i, j, k) = sum / num;
			}
			else if (den.isInBounds(Vec3i(i, j, k), 2) && flags.isSrc(i, j + 2, k) && vh(i, j + 2, k) != -2) {
				Real sum = 0;
				Real num = 0;
				for (int z = -1; z <= 1; z++) {
					for (int y = -1; y <= 1; y++) {
						for (int x = -1; x <= 1; x++) {
							if (vh(i + x, j + y + 1, k + z) >= 0) {
								sum += den(i + x, j + y + 1, k + z);
								num++;
							}
						}
					}
				}
				if (num>0) den(i, j, k) = sum / num;
			}
		}
	}

	inline void addRhoRegsSrc(std::vector<Triplet<Real>>& triplets, const FlagGrid& vhInflow, const FlagGrid& flags, const Vec3i& s, const int dim, const Real smooth, const Real kinetic, const Real sigma) {
		// add rho and regs
		Vec3i strides = Vec3i(1, s.x, s.x*s.y);
		Real mult;
		for (int k = 0; k < s.z; k++) {
			for (int j = 0; j < s.y; j++) {
				for (int i = 0; i < s.x; i++) {
					if (!flags.isSrc(i, j, k)) continue;
					int ijk = vhInflow.index(i, j, k);
					int srcIdx = vhInflow(ijk);
					if (srcIdx < 0) continue;

					// smoothness and kinetic energy regularizer if neighbour in visual hull
					mult = 0;
					for (int c = 0; c < dim; c++) {
						int idxLow = vhInflow(ijk - strides[c]);
						int idxUp = vhInflow(ijk + strides[c]);
						if (idxLow >= 0 && flags.isSrc(idxLow)) {
							triplets.push_back(Triplet<Real>(srcIdx, idxLow, -smooth));
							mult++;
						}
						if (idxUp >= 0 && flags.isSrc(idxUp)) {
							triplets.push_back(Triplet<Real>(srcIdx, idxUp, -smooth));
							mult++;
						}
					}
					triplets.push_back(Triplet<Real>(srcIdx, srcIdx, sigma + mult*smooth + kinetic));
				}
			}
		}
	}

	KERNEL(bnd=2) void setupVH(FlagGrid& vh, Grid<Real>& den, const Grid<Real>& densityTarget, const MACGrid& vel, const FlagGrid& flags, const Real weightsThresh, const Vec3i& s, const Vec3i& strides, const int *neighbours) {
		if (!flags.isSrc(i, j, k) && flags.isFluid(i,j,k) && fabs(densityTarget(i, j, k)) < 1e-2) {
			vh(i, j, k) = -2; // mark cell in R as outside vh

			// trace back
			Vec3 currentVel = vel.getCentered(i, j, k);
			Vec3d pos = Vec3d(i + 0.5, j + 0.5, k + 0.5) - Vec3d(currentVel.x, currentVel.y, currentVel.z);

			// interpolate between neighbouring cells
			Real weights[8];
			int idx = getWeights(weights, pos, s, strides.z); // lowest cell
			for (int ind = 0; ind < 8; ind++) {
				if (flags.isSrc(idx + neighbours[ind]) && fabs(weights[ind]) > weightsThresh) {
					den(idx + neighbours[ind]) = 0; // delete inflow that would lead to zero target cell
					vh(idx + neighbours[ind]) = -2; // mark cell in RI as outside vh if not already in the system
				}
			}
		}
	}

	KERNEL(single, bnd = 2) void setupIndex(FlagGrid& vh, Vec3i& Ns, const MACGrid& vel, const FlagGrid& flags, const Vec3i& s, const Vec3i& strides, const int *neighbours) {
		if (!flags.isSrc(i, j, k) && flags.isFluid(i, j, k) && vh(i, j, k) != -2) {
			// trace back
			Vec3 currentVel = vel.getCentered(i, j, k);
			Vec3d pos = Vec3d(i + 0.5, j + 0.5, k + 0.5) - Vec3d(currentVel.x, currentVel.y, currentVel.z);

			// interpolate between neighbouring cells
			bool oneValidSrcCell = false;
			Real weights[8];
			int idx = getWeights(weights, pos, s, strides.z); // lowest cell
			for (int ind = 0; ind < 8; ind++) {
				if (flags.isSrc(idx + neighbours[ind]) && vh(idx + neighbours[ind]) != -2 /*&& fabs(weights[ind])>1e-4*/) {
					oneValidSrcCell = true;
					if(vh(idx + neighbours[ind]) == -1) vh(idx + neighbours[ind]) = Ns.y++; // index cell in RI 
				}
			}
			if (oneValidSrcCell) vh(i, j, k) = Ns.x++; // index cell in R 
		}
	}

	KERNEL(single, bnd=2) void setupA(FlagGrid& vh, std::vector<Triplet<Real>>& triplets, const MACGrid& vel, const FlagGrid& flags, const Vec3i& s, const Vec3i& strides, const int *neighbours) {
		if (!flags.isSrc(i, j, k) && vh(i,j,k)>=0) {
			// trace back
			Vec3 currentVel = vel.getCentered(i, j, k);
			Vec3d pos = Vec3d(i + 0.5, j + 0.5, k + 0.5) - Vec3d(currentVel.x, currentVel.y, currentVel.z);

			// interpolate between neighboring cells
			Real weights[8];
			int idx = getWeights(weights, pos, s, strides.z); // lowest cell
			for (int ind = 0; ind < 8; ind++) {
				if (flags.isSrc(idx + neighbours[ind]) /*&& fabs(weights[ind])>1e-4*/ && vh(idx + neighbours[ind])>=0) {
					triplets.push_back(Eigen::Triplet<Real>(vh(i, j, k), vh(idx + neighbours[ind]), weights[ind]));
				}
			}
		}
	}

	KERNEL(bnd = 2, reduce = +) returns(Real sum = 0.0)
		Real knSumMissingDensity(const Grid<Real>& denResidual, const FlagGrid& vh, const FlagGrid& flags) {
		if (!flags.isSrc(i, j, k) && flags.isFluid(i, j, k) && vh(i, j, k) < 0) sum += denResidual(i, j, k);
	}

	KERNEL(bnd = 2, reduce = +) returns(Real sum = 0.0)
		Real knSumMatchedDensity(const Grid<Real>& denResidual, const FlagGrid& vh, const FlagGrid& flags) {
		if (!flags.isSrc(i, j, k) && flags.isFluid(i, j, k) && vh(i, j, k) >= 0) sum += denResidual(i, j, k);
	}

	KERNEL(bnd = 2, reduce = max) returns(Real absMaxNegValue = 0.0)
		Real knAbsMaxNegValue(const Grid<Real>& denResidual, const FlagGrid& vh, const FlagGrid& flags) {
		if (!flags.isSrc(i, j, k) && flags.isFluid(i, j, k) && vh(i, j, k) >= 0) {
			if (denResidual(i, j, k) < 0 && fabs(denResidual(i, j, k)) > absMaxNegValue) absMaxNegValue = fabs(denResidual(i, j, k));
		}
	}

	KERNEL(single, idx) void setupB(const Grid<Real>& denResidual, const Grid<Real>& denAdvect, const FlagGrid& vh, const FlagGrid& flags, const Real portionMissingDensity, const Real absMaxNegValue, VectorX& b) {
		if (!flags.isSrc(idx) && vh(idx) >= 0) {
			//b[vh(idx)] = max(denResidual(idx) + avgMissingDensity, -denAdvect(idx));
			b[vh(idx)] = max(denResidual(idx) + (denResidual(idx) + absMaxNegValue)*portionMissingDensity, -denAdvect(idx));
		}
	}

	KERNEL(idx) void ensureNNInflow(const FlagGrid& vh, const FlagGrid& flags, VectorX& z, const Grid<Real>& initialDensity) {
		if (vh(idx) >= 0 && flags.isSrc(idx)) {
			if (z[vh(idx)] + initialDensity(idx) >= 0) return;
			z[vh(idx)] = -initialDensity(idx);
		}
	}

	KERNEL(idx, reduce = +) returns(int sum = 0)
	int countSrcCells(const FlagGrid& flags, int a) {
		if(flags.isSrc(idx)) sum ++;
	}

	// estimate update of values in src in den in order to match density in imgs after advecting with vel
	void estimateDensityInflow(Grid<Real>& den, const MACGrid& vel, const FlagGrid& flags, TomoParams& tomoParams, const OFParams& ofParams, const Grid<Real>& densityTarget) {
		FluidSolver* solver = den.getParent();
		
		// create visual hull and A
		FlagGrid vh(solver); vh.setConst(-1);
		Grid<Real> densityTarRes(densityTarget);
		Vec3i Ns = Vec3i(0); // Ns.x = number of target cells (\Gamma_R), Ns.y = number of src cells (\Gamma_{RI})
		Vec3i strides = ofParams.getStrides();
		const int neighbours[8] = { 0 , strides.y, strides.x, strides.x + strides.y, strides.z, strides.y + strides.z, strides.x + strides.z, strides.x + strides.y + strides.z };
		setupVH(vh, den, densityTarRes, vel, flags, 1e-5, den.getSize(), strides, neighbours);
		setupIndex(vh, Ns, vel, flags, den.getSize(), strides, neighbours);
		delOutsideSrc(den, vh, flags); 

		if (Ns.x == 0) errMsg("Zero target cells in estimateDensityInflow!!!");
		if (Ns.y == 0) errMsg("Zero source cells in estimateDensityInflow!!!");
		std::vector<Triplet<Real>> triplets;
		triplets.reserve(8 * Ns.x);
		setupA(vh, triplets, vel, flags, den.getSize(), strides, neighbours);

		// create densityResidual
		Grid<Real> densityAdvect(den);
		advectSemiLagrange(&flags, &vel, &densityAdvect, ORDER, MACCORMACK_STRENGTH);
		densityTarRes.sub(densityAdvect);
		//Real avgMissingDensity = knSumMissingDensity(densityResidual, vh, flags) / Ns.x;
		Real absMaxNegValue = knAbsMaxNegValue(densityTarRes, vh, flags);
		Real portionMissingDensity = knSumMissingDensity(densityTarRes, vh, flags) / (knSumMatchedDensity(densityTarRes, vh, flags) + Ns.x*absMaxNegValue);

		int numSrcCells = countSrcCells(flags, 0);
		Real avgMissingDensity = (knSumMissingDensity(densityTarRes, vh, flags) + knSumMatchedDensity(densityTarRes, vh, flags)) / numSrcCells;
		//if(avgMissingDensity>0) setSrcToConst(den, flags, avgMissingDensity);
		//return;

		// setup LSE and solver
		LeastSquaresConjugateGradient<SparseMatrix<Real>> m_lscg;
		SparseMatrix<Real, RowMajor> m_P;
		SparseMatrix<Real, ColMajor> m_reg;

		m_P.resize(Ns.x, Ns.y);
		m_P.setFromTriplets(triplets.begin(), triplets.end());

		triplets.clear();
		triplets.reserve(5 * Ns.y);
		addRhoRegsSrc(triplets, vh, flags, den.getSize(), 3, tomoParams.getSmoothInflow(), tomoParams.getKineticInflow(), tomoParams.getSigma());
		m_reg.resize(Ns.y, Ns.y);
		m_reg.setFromTriplets(triplets.begin(), triplets.end());
		m_regLSCG = &m_reg;

		m_lscg.compute(m_P);
		if (m_lscg.info() != Success) errMsg("tomography(): compute P failed");
		m_lscg.setTolerance(1e-2);

		VectorX x; x.setZero(Ns.y);
		VectorX b; b.setZero(Ns.x);
		//setupB(densityResidual, densityAdvect, vh, flags, avgMissingDensity, b);
		setupB(densityTarRes, densityAdvect, vh, flags, portionMissingDensity, absMaxNegValue, b);

		m_lscg.setMaxIterations(1000 * max(Ns.x, Ns.y));

		VectorX z; z.setZero(Ns.y);
		VectorX y; y.setZero(Ns.y);
		VectorX zPrev; zPrev.setZero(Ns.y);
		VectorX m_rhsUpdate; m_rhsUpdate.setZero(Ns.y);
		for (int iter = 0; iter < 10; iter++) {
			// x-update
			m_rhsUpdate = x + tomoParams.getSigma()*y;
			m_rhsUpdateLSCG = &m_rhsUpdate;
			x += tomoParams.getSigma() * (y - m_lscg.solveWithGuess(b, z));

			// z-update
			z -= tomoParams.getTau()*x;
			ensureNNInflow(vh, flags, z, den);

			// y-update
			y = z + tomoParams.getTheta()*(z - zPrev);

			zPrev = z;
		}
		x = z;

		updateSrc(den, vh, flags, x);
		extrapolateSrc(den, vh, flags);

		m_regLSCG = nullptr;
	}

	PYTHON() void predict(MACGrid& vel, Grid<Real>& denPredict, const Grid<Real>& den, Grid<Real>& pressure, FlagGrid& flags, TomoParams* tomoParams, const OFParams& ofParams, const Grid<Real>* densityTarget, const Real cgAccuracy = 1e-3, const Real alpha = 0) {
		// advect velocity with itself
		advectSemiLagrange(&flags, &vel, &vel, ORDER, MACCORMACK_STRENGTH);
		solvePressure(vel, pressure, flags, cgAccuracy, 0, 0, 0, 1e-4, 1.5, true, 3, false, false, true);
		if (alpha > 0) {
			setWallBcs(flags, vel);
			cgSolveDiffusion(flags, vel, alpha);
		}

		// advect density forward
		denPredict.copyFrom(den);
		if (vel.getMaxAbs() > 0) {
			if (tomoParams) estimateDensityInflow(denPredict, vel, flags, *tomoParams, ofParams, *densityTarget); // estimate source
			advectSemiLagrange(&flags, &vel, &denPredict, ORDER, MACCORMACK_STRENGTH);
			deleteSmokeInBndsAdaptFlagsKN(denPredict, flags);
		}
	}
	
	PYTHON() void align(MACGrid& vel, Grid<Real>& den, Grid<Real>& pressure, const MACGrid& velUpdate, FlagGrid& flags, TomoParams* tomoParams, const OFParams& ofParams, const Grid<Real>* densityTarget, const Real cgAccuracy = 1e-3) {
		// add update
		vel.add(velUpdate);

		// estimate source
		if (tomoParams) estimateDensityInflow(den, vel, flags, *tomoParams, ofParams, *densityTarget);
		//setSrcToConst(src0, flags, 500);

		advectSemiLagrange(&flags, &vel, &den, ORDER, MACCORMACK_STRENGTH);
		deleteSmokeInBndsAdaptFlagsKN(den, flags);
	}

	// -------------------- down- and upsample for multi-scale approaches --------------------

	inline void findSrc(OFParams& ofParams, const FlagGrid& flagsC) {
		bool foundSrc = false;
		int centerX = flagsC.getSizeX() / 2;
		int centerZ = flagsC.getSizeZ() / 2;
		for (int j = 0; j <flagsC.getSizeY(); j++) {
			if (!foundSrc && flagsC.isSrc(centerX, j, centerZ)) foundSrc = true;
			if (foundSrc && !flagsC.isSrc(centerX, j, centerZ)) {
				ofParams.setHeightOfSrc(j - 1);
				break;
			}
		}
		if (ofParams.isBivariate() && !foundSrc) errMsg("Error in downsample, did not find src flag.");
	}

	inline void downsample(Grid<Real>& denPredictC, Grid<Real>& den1TC, Grid<Real>* denC, Grid<Real>* denMaskC, FlagGrid& flagsC, MACGrid* velPredictC, OFParams& ofParams, const Grid<Real>& denPredict, const Grid<Real>& den1T, const Grid<Real>* den, const Grid<Real>* denMask, const FlagGrid& flags, const MACGrid* velPredict) {
		// adapt params depending on grid size
		ofParams.updateSize(denPredictC.getSize());

		// interpolate flagGrid; make sure to set obstacle cells to zero in density grids
		downsampleFlagGrid(flagsC, flags); findSrc(ofParams, flagsC);

		// interpolate denPredict
		interpolateGrid(denPredictC, denPredict); deleteSmokeInBndsAdaptFlagsKN(denPredictC, flagsC);
		interpolateGrid(den1TC, den1T); deleteSmokeInBndsAdaptFlagsKN(den1TC, flagsC);
		
		if (denC) {
			interpolateGrid(*denC, *den); deleteSmokeInBndsAdaptFlagsKN(*denC, flagsC);
		}

		if (denMaskC) {
			interpolateGrid(*denMaskC, *denMask); deleteSmokeInBndsAdaptFlagsKN(*denMaskC, flagsC);
		}

		if (velPredict && velPredictC) {
			interpolateMACGrid(*velPredictC, *velPredict);
			velPredictC->multConst(0.5);
		}
	}

	KERNEL() void deleteSrcLowerOutflowCells(MACGrid& velUpdate, const FlagGrid& flags) {
		if (flags.isSrc(i,j,k) || j<2) velUpdate(i,j,k) = Vec3(0.);
	}

	inline void upsample(MACGrid& velUpdate, OFParams& ofParams, const MACGrid& velUpdateC, const FlagGrid& flags) {
		// adapt params depending on grid size
		ofParams.updateSize(velUpdate.getSize());

		// interpolate
		interpolateMACGrid(velUpdate, velUpdateC);
	
		// scale velocity
		velUpdate.multConst(2.);
	}

	// --------------------------------------------------------------------------------------
	// -------------------- optical flow part --------------------
	// --------------------------------------------------------------------------------------

	inline void calcVelUpdate(MACGrid& velUpdate, const MACGrid* velCoarse, const Grid<Real>& den0, const Grid<Real>& den1, FlagGrid& flags, Grid<Real>& pressure, OFParams& ofParams) {
		myClock::time_point startTime = myClock::now();

		DivFreeOF of(velUpdate.getParent(), ofParams, flags, velCoarse, den0, pressure, &den1, nullptr);
		of.solve(velUpdate, velCoarse, flags);

		Real denError;
		Real eqError = DivFreeOFNS::computeError(denError, flags, velUpdate, den0, ofParams, &den1);
		//debMsg("t=" << ofParams.getT() << ", res=" << flags.getSizeX() << ": calcVelUpdate() took " << diffClock(myClock::now(), startTime) << " ms, max vel=" << velUpdate.getMaxAbs() << ", eq err=" << eqError << ", adv err=" << denError, 1);
	}

	void calcVelUpdateMS(MACGrid& velUpdate, const Grid<Real>& den0, const Grid<Real>& den1, FlagGrid& flags, Grid<Real>& pressure, OFParams& ofParams) {
		MACGrid* velUpdateLowerUpsampled = nullptr;

		if (velUpdate.getSizeX() > ofParams.getMinSize()) {
			velUpdateLowerUpsampled = new MACGrid((velUpdate.getParent()));

			// downsample
			ofParams.scaleSmooth(0.5); ofParams.scaleKinetic(0.5);

			FluidSolver solverLower = FluidSolver(velUpdate.getSize() / 2, ofParams.getDim(), -1);
			MACGrid velUpdateLower(&solverLower); // stays empty, is the return value of the recursive function call, no downsample, only upsample
			Grid<Real> den0Lower(&solverLower); Grid<Real> den1Lower(&solverLower); Grid<Real> pressureLower(&solverLower); FlagGrid flagsLower(&solverLower);
			downsample(den0Lower, den1Lower, nullptr, nullptr, flagsLower, nullptr, ofParams, den0, den1, nullptr, nullptr, flags, nullptr);
			// recursive call
			calcVelUpdateMS(velUpdateLower, den0Lower, den1Lower, flagsLower, pressureLower, ofParams);

			// upsample
			ofParams.scaleSmooth(2.); ofParams.scaleKinetic(2.);
			upsample(*velUpdateLowerUpsampled, ofParams, velUpdateLower, flags);
			// not needed in univariate case
			//solvePressure(*velUpdateLowerUpsampled, pressure, flags, ofParams.getPressureAcc());
		}

		// move den0 forward with coarse velocity estimate
		if (velUpdateLowerUpsampled) {
			Grid<Real> denPredictNew(den0);
			advectSemiLagrange(&flags, velUpdateLowerUpsampled, &denPredictNew, ORDER, MACCORMACK_STRENGTH);
			deleteSmokeInBndsAdaptFlagsKN(denPredictNew, flags);
			calcVelUpdate(velUpdate, velUpdateLowerUpsampled, denPredictNew, den1, flags, pressure, ofParams);
			// combine velUpdate from coarser scales with current one
			velUpdate.add(*velUpdateLowerUpsampled);
			delete velUpdateLowerUpsampled;
		}
		else {
			calcVelUpdate(velUpdate, velUpdateLowerUpsampled, den0, den1, flags, pressure, ofParams);
		}
	}

	// optical flow on staggered grid
	// includes obstacle and bnd handling
	//		obstacles do not move, no fractions supported (yet), no empty cells supported, only slip bcs
	PYTHON() void reconstructVelocity(MACGrid& velUpdate, const Grid<Real>& den0, const Grid<Real>& den1, FlagGrid& flags, Grid<Real>& pressure, OFParams& ofParams) {
		myClock::time_point startTime = myClock::now();

		velUpdate.clear();
		calcVelUpdateMS(velUpdate, den0, den1, flags, pressure, ofParams);

		// final error calculation
		Real denError;
		Real eqError = DivFreeOFNS::computeError(denError, flags, velUpdate, den0, ofParams, &den1);
		//debMsg("t=" << ofParams.getT() << ": whole reconVel() took " << diffClock(myClock::now(), startTime) << " ms, eq err=" << eqError << ", den err=" << denError, 1);
	 }


	// --------------------------------------------------------------------------------------
	// -------------------- tomography part --------------------
	// --------------------------------------------------------------------------------------

	/* not imported from old code:
	*  - delete upper part of density in volume to simulate outflow
	*  - interpolate image grid
	*  - regular perspective cameras with focal length etc. (pvFactor, distVolImg), now: can't rotate around arbitrary angle, only 3 orthogonal cameras or real cameras with ray files
	*  - solve for source: setup matrix system that has only one row, used to adjust overall brightness to sum up to match denAmountImg,
	*    some form of scaling; done with cgls, mostly important for single scale (i.e. where density in smoke area expands in depth but new source gets added)
	*/

	// calculate density based on images in imgs and ray information in i, consider visual hull
	// if solving for density update correction, imgs is imgs-imgsPredicted-imgsDenUpate
	PYTHON() void reconstructDensity(Grid<Real>& density, const Image& i, const Grid<Real>& imgs, const FlagGrid& flags, TomoParams& tomoParams, const Grid<Real>* densityMask = nullptr) {
		myClock::time_point startTime = myClock::now();

		Tomography tomo(tomoParams, imgs, i, flags, densityMask);
		tomo.solve(density, imgs, nullptr, nullptr);

		debMsg("t="<< tomoParams.getT()<<", reconDensity() took " << diffClock(myClock::now(), startTime) << " ms, max den=" << density.getMaxAbs() << ", tomo error=" << tomo.computeError(density, imgs), 1);
	}

	// --------------------------------------------------------------------------------------
	// -------------------- combined OF and tomography --------------------
	// --------------------------------------------------------------------------------------	

	PYTHON() int reconstructDenVel(MACGrid& velUpdate, const MACGrid& velPredict, const Image& i, const Grid<Real>& imgs, const Grid<Real>& denPredict, FlagGrid& flags, OFParams& ofParams, TomoParams& tomoParams, Grid<Real>& pressure, const Grid<Real>* densityMask=nullptr, const Grid<Real>* densityTarget=nullptr) {
		myClock::time_point startTime0 = myClock::now();

		Tomography tomo(tomoParams, imgs, i, flags, densityMask);
		printf("t=%03d,  tomo() took %07.f ms \n", ofParams.getT(), diffClock(myClock::now(), startTime0));

		DivFreeOF of(velUpdate.getParent(), ofParams, flags, nullptr, denPredict, pressure, densityTarget, &tomo);
		printf("t=%03d,    of() took %07.f ms \n", ofParams.getT(), diffClock(myClock::now(), startTime0));
		fflush(stdout);

		int iter = of.solve(velUpdate, &velPredict, flags, &denPredict, &imgs, &tomo, densityTarget);

		// final error calculation
		of.m_helper->add(denPredict); // denUpdate not needed anymore, used for error calculation
		Real denError;
		Real eqError = DivFreeOFNS::computeError(denError, flags, velUpdate, denPredict, ofParams, of.m_helper);
		printf("Single scale t=%03d, iter=%03d, eq err=%.4f , adv err=%.4f , tomo error=%.4f\n", ofParams.getT(), iter, eqError, denError, tomo.computeError(*of.m_helper, imgs));
		fflush(stdout); 
		return iter;
	}

	int recursiveCalcUpdate(MACGrid& velUpdate, const MACGrid& velPredict, const Image& i, const Grid<Real>& imgs, const Grid<Real>& denPredict, const Grid<Real>& den, FlagGrid& flags, OFParams& ofParams, TomoParams& tomoParams, Grid<Real>& pressure, const Grid<Real>& densityTarget, const Grid<Real>* densityMask) {
		int iter = 0;
		MACGrid* velUpdateLowerUpsampled = nullptr; 
		// go one level down
		if (denPredict.getSizeX() > ofParams.getMinSize()) {
			velUpdateLowerUpsampled = new MACGrid((velUpdate.getParent()));

			// downscale denP, densityMask, flags, adapt params
			// parameter adjustment
			int heightOfSrcUpper = ofParams.getHeightOfSrc(); 
			ofParams.scaleSmooth(0.5); ofParams.scaleKinetic(0.5); tomoParams.scaleSmooth(0.5); tomoParams.scaleKinetic(0.5);

			FluidSolver solverLower = FluidSolver(flags.getSize() / 2, ofParams.getDim(), -1);
			MACGrid velUpdateLower(&solverLower); Grid<Real> pressureLower(&solverLower);
			MACGrid velPredictLower(&solverLower); FlagGrid flagsLower(&solverLower); Grid<Real> denPredictLower(&solverLower); Grid<Real> denLower(&solverLower); Grid<Real> densityTargetLower(&solverLower); Grid<Real>* densityMaskLower = nullptr;
			if(densityMask) densityMaskLower = new Grid<Real>(&solverLower);
			downsample(denPredictLower, densityTargetLower, &denLower, densityMaskLower, flagsLower, &velPredictLower, ofParams, denPredict, densityTarget, &den, densityMask, flags, &velPredict);
			
			// recursive call: calculate velocity update on lower scale
			iter += recursiveCalcUpdate(velUpdateLower, velPredictLower, i, imgs, denPredictLower, denLower, flagsLower, ofParams, tomoParams, pressureLower, densityTargetLower, densityMaskLower);

			// upsample
			ofParams.setHeightOfSrc(heightOfSrcUpper); 
			ofParams.scaleSmooth(2.); ofParams.scaleKinetic(2.); tomoParams.scaleSmooth(2.); tomoParams.scaleKinetic(2.);
			
			upsample(*velUpdateLowerUpsampled, ofParams, velUpdateLower, flags);
			solvePressure(*velUpdateLowerUpsampled, pressure, flags, ofParams.getPressureAcc());

			if (densityMask) delete densityMaskLower;
		}

		// velUpdate is u2, then u2,Final
		if (velUpdateLowerUpsampled) {
			MACGrid velPredictNew(velPredict);
			velPredictNew.add(*velUpdateLowerUpsampled);
			Grid<Real> denPredictNew(den);
			estimateDensityInflow(denPredictNew, velPredictNew, flags, tomoParams, ofParams, densityTarget); // estimate density inflow

			advectSemiLagrange(&flags, &velPredictNew, &denPredictNew, ORDER, MACCORMACK_STRENGTH);
			deleteSmokeInBndsAdaptFlagsKN(denPredictNew, flags);

			iter += reconstructDenVel(velUpdate, velPredictNew, i, imgs, denPredictNew, flags, ofParams, tomoParams, pressure, densityMask, &densityTarget);

			// combine velUpdate from coarser scales with current one
			velUpdate.add(*velUpdateLowerUpsampled);
			delete velUpdateLowerUpsampled;
		}
		else {
			iter += reconstructDenVel(velUpdate, velPredict, i, imgs, denPredict, flags, ofParams, tomoParams, pressure, densityMask, &densityTarget);
		}
		return iter;
	}

	PYTHON() int reconstructDenVelMS(MACGrid& velUpdate, const MACGrid& velPredict, const Image& i, const Grid<Real>& imgs, const Grid<Real>& denPredict, const Grid<Real>& den, FlagGrid& flags, OFParams& ofParams, TomoParams& tomoParams, Grid<Real>& pressure, const Grid<Real>& densityTarget, const Grid<Real>* densityMask = nullptr) {
		if (tomoParams.shapeLimitspecified()) errMsg("Error in reconstructDenVelMS: shapeLimit given but downscaling isn't implemented for it yet.");
		
		velUpdate.clear();
		int iter = recursiveCalcUpdate(velUpdate, velPredict, i, imgs, denPredict, den, flags, ofParams, tomoParams, pressure, densityTarget, densityMask);

		Real denError;
		Real eqError = DivFreeOFNS::computeError(denError, flags, velUpdate, denPredict, ofParams, 0);
		printf("All scales t=%03d, iter=%03d, eq err=%.4f , adv err=%.4f \n", ofParams.getT(), iter, eqError, denError);
		fflush(stdout);
		return iter;
	}

	KERNEL() void copyPartOfGridKN(MACGrid& vel, const MACGrid& velO) {
		vel(i, j, k) = velO(i, j, k);
	}

	PYTHON() void copyPartOfGrid(MACGrid& vel, const MACGrid& velO) {
		copyPartOfGridKN(vel, velO);
	}

	KERNEL() void copyPartOfGridDenKN(Grid<Real>& den, const Grid<Real>& denO) {
		den(i, j, k) = denO(i, j, k);
	}
	PYTHON() void copyPartOfGridDen(Grid<Real>& den, const Grid<Real>& denO) {
		copyPartOfGridDenKN(den, denO);
	}

	KERNEL(idx, reduce = max) returns(Real maxVal = 0)
	Real getMaxYVel(const MACGrid& vel) {
		if (vel[idx].y > maxVal) maxVal = vel[idx].y;
	}

	KERNEL(idx, reduce = +) returns(Real sum = 0)
	Real getAvgYVelMasked(const MACGrid& vel, const Grid<Real>& den) {
		if (den(idx)>1e-2) sum += vel[idx].y;
	}

	KERNEL(idx, reduce = +) returns(Real sum = 0)
	Real getYVelMask(const Grid<Real>& den) {
		if (den(idx)>1e-2) sum ++;
	}

	PYTHON() void estimatePlumeSpeed(MACGrid& vel0, MACGrid& vel1, Grid<Real>& den, Grid<Real>& denPast, Grid<Real>& denFuture, const Image& i, const Grid<Real>& imgs, const Grid<Real>& imgsPast, const Grid<Real>& imgsFuture, FlagGrid& flags, Grid<Real>& pressure, TomoParams& tomoParams, OFParams& ofParams, const bool tomoExists) {
		if (tomoExists) {
			
		}
		else {
			Tomography tomo(tomoParams, imgs, i, flags, nullptr);
			tomo.solve(den, imgs, nullptr, nullptr);
		}

		if (tomoExists) {

		}
		else {
			Tomography tomoPast(tomoParams, imgsPast, i, flags, nullptr);
			tomoPast.solve(denPast, imgsPast, nullptr, nullptr);
		}

		if (tomoExists) {

		}
		else {
			Tomography tomoFuture(tomoParams, imgsFuture, i, flags, nullptr);
			tomoFuture.solve(denFuture, imgsFuture, nullptr, nullptr);
		}

		calcVelUpdateMS(vel0, denPast, den, flags, pressure, ofParams);
		calcVelUpdateMS(vel1, den, denFuture, flags, pressure, ofParams);

		
		denPast.add(den);
		denPast.multConst(0.5);
		denFuture.add(den);
		denFuture.multConst(0.5);

		if (getYVelMask(denPast) == 0) errMsg("empty denPast field!");
		if (getYVelMask(denFuture) == 0) errMsg("empty denFuture field!");
		debMsg(ofParams.getT()<<": max0, max1, avg0, avg1: "<< getMaxYVel(vel0) << ", " << getMaxYVel(vel1) << ", " << getAvgYVelMasked(vel0, denPast) / getYVelMask(denPast) << ", " << getAvgYVelMasked(vel1, denFuture) / getYVelMask(denFuture),1);

	}


} // namespace
