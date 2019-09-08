/******************************************************************************
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * @author: Marie-Lena Eckert http://marielenaeckert.com/
 *
 * data structures and functions useful to both tomography and divFreeOF PD solves
 *
******************************************************************************/

#ifndef PD_SOLVE
#define PD_SOLVE

#include <chrono>
#include "manta.h"
#include "commonkernels.h"
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>

// -------------------- global variables for LSCG solver --------------------

extern Eigen::SparseMatrix<Real, Eigen::ColMajor>* m_regLSCG;
extern Eigen::Matrix<Real, -1, 1>* m_rhsUpdateLSCG;

// -------------------- end global variables for LSCG solver --------------------


using myClock = std::chrono::high_resolution_clock;
using namespace std;

namespace Manta {

	// -------------------- constants and typedefs --------------------
	const Real MACCORMACK_STRENGTH = 0.8; 
	const int ORDER = 2;
	typedef Eigen::Matrix<Real, -1, 1> VectorX;

	// -------------------- externally defined functions --------------------
	extern void solvePressure(MACGrid& vel, Grid<Real>& pressure, const FlagGrid& flags, Real cgAccuracy = 1e-3,
		const Grid<Real>* phi = 0,
		const Grid<Real>* perCellCorr = 0,
		const MACGrid* fractions = 0,
		Real gfClamp = 1e-04,
		Real cgMaxIterFac = 1.5,
		bool precondition = true, // Deprecated, use preconditioner instead
		int preconditioner = 1,
		bool enforceCompatibility = false,
		bool useL2Norm = false,
		bool zeroPressureFixing = false,
		const Grid<Real> *curv = NULL,
		const Real surfTens = 0.,
		Grid<Real>* retRhs = NULL);

	extern void advectSemiLagrange(const FlagGrid* flags, const MACGrid* vel, GridBase* grid,
		int order = 1, Real strength = 1.0, int orderSpace = 1, int clampMode = 2);

	extern void projectPpmFull(const Grid<Real>& val, std::string name, int shadeMode = 0, Real scale = 1.);

	// -------------------- useful class definitions --------------------

	PYTHON()
		class ShapeDetails : public PbClass {
		public:
			PYTHON() ShapeDetails(FluidSolver *parent, const int shape, const Vec3 center, const Vec3 vec, const Real radius);
			const int shape; // 0:cylinder, 1:sphere, 2:box
			const Vec3 center;
			const Vec3 vec;
			const Real radius;
			
			PYTHON() int getHeightOfSrc() const;
	};

	PYTHON()
		class PDParams : public PbClass {
		public:
			PYTHON() PDParams(FluidSolver *parent, const Real sigma, const Real tau, const Real theta, const int mxIter);
			const Real sigma;
			const Real tau;
			const Real theta;
			const int maxIter;
	};

	PYTHON()
		class RegWeightsTomo : public PbClass {
		public:
			PYTHON() RegWeightsTomo(FluidSolver *parent, const Real smooth, const Real kinetic, const Real smoothInflow=0, const Real kineticInflow=0);

			void scaleSmooth(Real factor);
			void scaleKinetic(Real factor); 

			Real getSmooth() const;
			Real getKinetic() const;
			Real getSmoothInflow() const;
			Real getKineticInflow() const;

		private:
			Real smooth;
			Real kinetic;
			Real smoothInflow;
			Real kineticInflow;
	};

	PYTHON()
		class TomoParams : public PbClass {
		public:
			PYTHON() TomoParams(FluidSolver *parent, const int t, const string path, const Real threshVH, const Real threshMask, const Real stepSize, const int minNumCams, const int numPixels, const int numVoxels, const Real angleWeight, const PDParams& pdParams, const RegWeightsTomo& regWeights, const ShapeDetails* shapeLimit = nullptr);
			
			PYTHON() void setT(int newT);
			void setNumPixels(int newNumPixels);
			void setNumVoxels(int newNumVoxels);
			void scaleSmooth(Real factor);
			void scaleKinetic(Real factor);

			PYTHON() Real getSmooth() const;
			PYTHON() Real getKinetic() const;
			PYTHON() Real getSmoothInflow() const;
			PYTHON() Real getKineticInflow() const;
			PYTHON() Real getThreshVH() const;
			PYTHON() Real getThreshMask() const;
			int getT() const;
			string getPath() const;
			bool shapeLimitspecified() const;
			int getNumPixels() const;
			int getNumVoxels() const;
			bool zeroVoxels() const;
			Real getStepSize() const;
			Real getAngleWeight() const;
			Real getSigma() const;
			Real getTau() const;
			Real getTheta() const;
			int getMaxIter() const;
			int getMinNumCams() const;
			const ShapeDetails* getShapeLimit() const;

		private:
			int t;
			const string path;
			const Real threshVH;
			const Real threshMask;
			const Real stepSize;
			const int minNumCams;
			int numPixels;
			int numVoxels;
			const Real angleWeight;
			const PDParams pdParams;
			RegWeightsTomo regWeights;
			const ShapeDetails* shapeLimit;
			void printParams() const;
	};

	PYTHON()
		class RegWeightsOF : public PbClass {
		public:
			PYTHON() RegWeightsOF(FluidSolver *parent, const Real smooth, const Real kinetic, const Real kineticZ = -1, const bool adaptiveKineticZ = false, const Real sumZero = 0);

			void scaleSmooth(Real factor);
			void scaleKinetic(Real factor);

			Real getSmooth() const;
			Real getKinetic() const;
			Real getKineticZ() const;
			bool useAdaptiveKineticZ() const;
			Real getSumZero() const;

		private:
			Real smooth;
			Real kinetic;
			Real kineticZ;
			bool adaptiveKineticZ;
			const Real sumZero;
	};

	PYTHON()
		class OFParams : public PbClass {
		public:
			PYTHON() OFParams(FluidSolver *parent, const int t, const int minSize, const bool bivariate, const string path, const Vec3i s, const Vec3i strides, const int N, const int dim, const int heightOfSrc, const bool useDenTarget, const Real inflowValue, const PDParams& pdParams, const RegWeightsOF& regWeights);
			
			PYTHON() void setT(int newT);
			void setPressureAcc(Real newPressureAcc);
			void scalePressureAcc(Real factor);
			void updateSize(Vec3i newSize);
			void setHeightOfSrc(int newHeightOfSrc);
			void setNumOfVelCmp(int newNumOfVelCmp);
			void scaleSmooth(Real factor);
			void scaleKinetic(Real factor);

			PYTHON() Real getSmooth() const;
			PYTHON() Real getKinetic() const;
			int getT() const;
			int getMinSize() const;
			bool isBivariate() const;
			string getPath() const;
			Real getPressureAcc() const;
			int getS(int c) const;
			Vec3i getStrides() const;
			int getStride(int c) const;
			int getDim() const;
			int getN() const;
			int getHeightOfSrc() const;
			int getNumOfVelCmp() const;
			bool shouldUseDenTarget() const;
			Real getSigma() const;
			Real getTau() const;
			Real getTheta() const;
			int getMaxIter() const;
			Real getSumZero() const;
			Real getKineticZ() const;
			bool useAdaptiveKinZ() const;
			PYTHON() Real getInflowValue() const;
			
		private:
			int t;
			const int minSize;
			const bool bivariate;
			const string path;
			Real pressureAcc;
			Vec3i s;
			Vec3i strides;
			int N;
			const int dim;
			int heightOfSrc;
			int numOfVelCmp;
			bool useDenTarget;
			Real inflowValue;
			const PDParams pdParams;
			RegWeightsOF regWeights;
			void printParams() const;
	};


	// -------------------- general helper functions --------------------

	inline float diffClock(myClock::time_point clock1, myClock::time_point& clock2) {
		float time = (std::chrono::duration_cast<std::chrono::milliseconds>(clock1 - clock2)).count();
		clock2 = clock1;
		return time;
	}

	// get max component of error vector
	KERNEL(idx, reduce = max) returns(Real maxVal = 0.)
	Real getMaxError(const Grid<Real>& density, const VectorX& error) {
		if (error[idx] > maxVal) maxVal = error[idx];
	}

	// get sum of error vector
	KERNEL(idx, reduce = +) returns(Real sum = 0.0)
	Real getSumError(const Grid<Real>& density, const VectorX& error) {
		sum += error[idx];
	}

	// get max component of vector
	KERNEL(idx, reduce = max) returns(Real maxVal = 0.)
	Real getMaxVec(const FlagGrid& vh, const VectorX& vec, const int N) {
		if (vh(idx) >= 0 && fabs(vec[N+vh(idx)]) > maxVal) maxVal = fabs(vec[N+vh(idx)]);
	}

	// get sum of vector
	KERNEL(idx, reduce = +) returns(Real sum = 0.0)
	Real getSumVec(const FlagGrid& vh, const VectorX& vec, const int N) {
		if (vh(idx) >= 0) sum += vec[N + vh(idx)];
	}

}

#endif