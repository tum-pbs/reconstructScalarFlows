/******************************************************************************
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * @author: Marie-Lena Eckert http://marielenaeckert.com/
 *
 * class for tomography solver
 *
******************************************************************************/

#include "tomography.h"

using namespace std;

namespace Manta {

	int Tomography::getVHEntry(int idx) const {
		return (*m_vh)[idx];
	}

	Real Tomography::computeError(const Grid<Real>& density, const Grid<Real>& imgs) {
		// error from equation
		int N = imgs.getSizeX()*imgs.getSizeY()*imgs.getSizeZ();
		VectorX error; error.setZero(N);
		VectorX densityVec(m_params->getNumVoxels());
		TomographyNS::copyGridToVec(*m_vh, densityVec, density);
#if USECG
		TomographyNS::calcError(m_vhP, error, m_P_single*densityVec, imgs);
#else
		TomographyNS::calcError(m_vhP, error, m_P*densityVec, imgs);
#endif
		Real retError = getMaxError(imgs, error);
		//debMsg("Tomo Error:, max=" << retError << ", avg=" << getSumError(imgs, error) / N, 1);
		return retError;
	}

	Tomography::~Tomography() {
		delete m_vh;
		delete m_vhP;
	}

	Tomography::Tomography(TomoParams& params, const Grid<Real>& imgs, const Image& i, const FlagGrid& flags, const Grid<Real>* densityMask) {
		// init visual hull grids
		m_vh = new FlagGrid(flags.getParent());
		m_vhP = new FlagGrid(imgs.getParent());
		m_params = &params;
		m_i = &i;
		// init matrices, b, rhs, pd variables, compute solvers; densityBlurred only used for mask in vh
		init(flags, imgs, densityMask);
	}

	// init before each new solve
	inline void Tomography::init(const FlagGrid& flags, const Grid<Real>& imgs, const Grid<Real>* densityMask) {
		// init visual hull grids
		// fill P, PP+rhoI and set vh, vhP, numPixels, numVoxels; imgs is only used for visual hull here
		m_i->setupMatrix(*m_vh, *m_vhP, m_P, m_P_single, m_reg, flags, imgs, *m_params, USECG, densityMask);
		debMsg("Tomo setupMatrix(): #pixels in vh: " << m_params->getNumPixels() << ", #voxels in vh: " << m_params->getNumVoxels(), 1);
		if (m_params->zeroVoxels()) errMsg("Error in tomography: zero voxels in visual hull.");
		
		// init vectors
		m_b.setZero(m_params->getNumPixels());
		m_rhsUpdate.setZero(m_params->getNumVoxels());
		m_z.setZero(m_params->getNumVoxels());

		// compute matrix; init solver
		m_regLSCG = &m_reg;
		m_lscg.compute(m_P);
		if (m_lscg.info() != Success) errMsg("tomography(): compute P failed");
		m_lscg.setTolerance(1e-2);
	}

	// imgs is passed by value, is modified but should not affect caller method
	void Tomography::solve(Grid<Real>& density, const Grid<Real>& imgs, const Grid<Real>* initialDen, const FlagGrid* flags) {
		myClock::time_point startTime0 = myClock::now(); 
		m_regLSCG = &m_reg;
		
		// init vectors
		m_x.setZero(m_params->getNumVoxels());
		m_zPrev = m_z;
		m_y = m_z;

#if USECG
		m_b.setZero(m_params->numPixels);
#endif

		// need to modify right hand side if already initial density values (when solving for a density update)
		if (initialDen) {
			Grid<Real> imgsCorrect(imgs);
			m_i->render(imgsCorrect, *initialDen, m_params->getStepSize(), "", "", true, flags);
			// fill b with imgs
			TomographyNS::setupB(*m_vhP, m_b, imgsCorrect, m_params->getAngleWeight());
		}
		else {
			// fill b with imgs
			TomographyNS::setupB(*m_vhP, m_b, imgs, m_params->getAngleWeight());
		}

#if USECG
		m_b = m_P_single.transpose()*m_b;
#endif

		// start PD solve
		Real mxCut = 0;
		Real mxCutP = 0;
		int iter, numOfCGIter;
		myClock::time_point startTime;
		Real times[] = { 0.,0.,0. }; 
		for (iter = 0; iter < m_params->getMaxIter(); iter++) {
			// solve tomography
			startTime = myClock::now();
			// x-update
			m_rhsUpdate     = m_x + m_params->getSigma()*m_y;
			m_rhsUpdateLSCG = &m_rhsUpdate;
			m_x += m_params->getSigma() * (m_y - m_lscg.solveWithGuess(m_b, m_z));
			numOfCGIter = m_lscg.iterations();
			times[0] += diffClock(myClock::now(), startTime);
				
			// z-update
			m_z -= m_params->getTau()*m_x;
			mxCut = TomographyNS::ensureNN(*m_vh, m_z, initialDen);
			times[1] += diffClock(myClock::now(), startTime);
				
			// y-update
			m_y = m_z + m_params->getTheta()*(m_z - m_zPrev);
			
			// check stopping criterion
			m_zPrev = m_z - m_zPrev;
			Real thresh = 1e-5*getMaxVec(*m_vh, m_z, 0);
			if (mxCut < thresh || std::fabs(mxCut - mxCutP) < 1e-2*thresh) break;

			m_zPrev = m_z;
			mxCutP = mxCut;
		}

		TomographyNS::copyVecToGridZero(density, *m_vh, m_z);
		m_regLSCG = nullptr;
		m_rhsUpdateLSCG = nullptr;
		//printf("tomo i=%03d, denMx=%.3f, denAvg=%.3f, mxCut=%.4f, solve() took %07.f ms, tomo: %07.f ms, nn: %07.f ms, blur: %07.f ms \n", iter, density.getMaxAbs(), getSumVec(*m_vh, m_z) / m_params->numVoxels, mxCut, diffClock(myClock::now(), startTime0), times[0], times[1], times[2]);
	}
}