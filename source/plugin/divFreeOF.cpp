/******************************************************************************
* MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * @author: Marie-Lena Eckert http://marielenaeckert.com/
 *
 * class for divergence-free optical flow solver
 *
******************************************************************************/

#include "divFreeOF.h"

using namespace std;

namespace Manta {

	void DivFreeOF::setupSystem(const MACGrid* velCoarse, const FlagGrid& flags, const Grid<Real>& den, const Grid<Real>* den1T, const Tomography* tomo) {
		// reserve triplets
		int numEntries = m_params->getDim() * m_params->getDim() * m_params->getNumOfVelCmp();
		if (m_params->isBivariate()) numEntries += (m_params->getDim() * m_params->getDim() + 1) * m_numOfDenComp + /*4 * */ m_params->getDim() * m_numOfDenComp;
		vector<Triplet<Real>> triplets;
		triplets.reserve(numEntries);

		// fill triplets and b; fill m_A with a large tripletsvector takes most time
		DivFreeOFNS::setupA(*m_vh, triplets, flags, den, den1T, *m_params, m_vhDen, tomo);
		m_A.setFromTriplets(triplets.begin(), triplets.end());

		if (!m_params->isBivariate()) DivFreeOFNS::setupB(*m_vh, m_b, velCoarse, flags, den, *m_params, *den1T);
	}

	inline void DivFreeOF::verifyArguments(const Grid<Real>* den1, const Tomography* tomo) {
		if (m_params->isBivariate()) {
			if (!tomo) errMsg("Error in OF(): bivariate case but no tomography object is passed.");
			if (m_params->shouldUseDenTarget() && !den1) errMsg("Error in OF(): bivariate case, should use denTarget, but no den1T passed.");
		}
		else {
			if (!den1) errMsg("Error in OF(): univariate case but no grid den1 is passed.");
		}
	}

	DivFreeOF::~DivFreeOF() {
		delete m_vh;
		if (m_params->isBivariate()) delete m_vhDen;
	}

	DivFreeOF::DivFreeOF(FluidSolver* s, OFParams& params, const FlagGrid& flags, const MACGrid* velCoarse, const Grid<Real>& den, Grid<Real>& pressure, const Grid<Real>* den1, const Tomography* tomo) {
		m_params = &params;

		// check if correct arguments for univariate and bivariate case
		verifyArguments(den1, tomo);

		init(s, flags, velCoarse, den, pressure, den1, tomo);
	}

	// takes its time (res=60, 500ms)
	void DivFreeOF::init(FluidSolver* s, const FlagGrid& flags, const MACGrid* velCoarse, const Grid<Real>& den, Grid<Real>& pressure, const Grid<Real>* den1T, const Tomography* tomo) {
		m_helper = &pressure;
		m_vh = new Grid<Vec3i>(s);
		// setup vh for density (larger than vh for tomography)
		m_numOfDenComp = 0;
		if (m_params->isBivariate()) {
			m_vhDen = new FlagGrid(s);
			m_numOfDenComp = DivFreeOFNS::setupDenVH(*m_vhDen, flags, den, tomo);
		}

		// setup visual hull for velocities
		m_params->setNumOfVelCmp(DivFreeOFNS::setupVelVH(*m_vh, flags, m_params->getDim()));
		m_size = m_params->getNumOfVelCmp() + m_numOfDenComp;
		debMsg("DivFreeOF: #vel in vh: " << m_params->getNumOfVelCmp() << ", #den in vh: " << m_numOfDenComp, 1);
		if (m_params->isBivariate() && m_numOfDenComp == 0) errMsg("Error in combined reconstruction: zero voxels in density visual hull.");

		// init vectors and matrix
		m_b.setZero(m_size);
		m_x.setZero(m_size);
		m_z.setZero(m_size);
		m_zPrev.setZero(m_size);
		m_y.setZero(m_size);
		m_A.resize(m_size, m_size);

		m_params->setPressureAcc(1e-2);

		// fill A, takes its time (res=60, 300ms)
		setupSystem(velCoarse, flags, den, den1T, tomo);

		// compute matrix; init solver
		m_cg.compute(m_A);
		if (m_cg.info() != Success) errMsg("divFreeOF(): compute A failed");
		m_cg.setTolerance(1e-2);
		m_cg.setMaxIterations(3e-3*m_cg.maxIterations());
	}

	KERNEL() void copyVelIntoSrc(MACGrid& velUpdate, const FlagGrid& flags, const int heightOfSrc) {
		if (flags.isSrc(i, j, k)) velUpdate(i, j, k).y = velUpdate(i, heightOfSrc + 1, k).y;
	}

	KERNEL(bnd=2) void setVelInflowAtFace(MACGrid& vel, FlagGrid& flags, const Real value, const Grid<Real>* denPredict, const Grid<Real>& denUpdate) {
		if (!flags.isSrc(i,j,k) && flags.isSrc(i, j - 1, k)) {
			Real velValue = max(vel(i, j, k).y, value);
			Vec3 pos = Vec3(i + 0.5, j + velValue, k + 0.5);
			Real denValue = denPredict ? denPredict->getInterpolated(pos) + denUpdate.getInterpolated(pos) : 1;
			if (fabs(denValue) > 1e-3 || flags.isInflow(i, j, k)) {
				vel(i, j, k).y = velValue;
				flags(i, j, k) |= FlagGrid::TypeInflow;
			}
		}
		else if (flags.isSrc(i, j, k)) vel(i, j, k) = Vec3(0.);
	}

	KERNEL(idx) void setSrcToObs(FlagGrid& flagsPressure, int a) {
		if (flagsPressure.isSrc(idx)) flagsPressure(idx) = FlagGrid::TypeObstacle;
	}

	void DivFreeOF::makeDivFree(MACGrid& velUpdate, const MACGrid* velPredict, FlagGrid& flags, const FlagGrid& flagsPressure, const Grid<Real>* denPredict, const Real velInflowValue, const bool setInflow, const bool bivariate) {
		if (velPredict) velUpdate.add(*velPredict);
		if (setInflow) {
			if (bivariate) setVelInflowAtFace(velUpdate, flags, velInflowValue, denPredict, *m_helper);
			solvePressure(velUpdate, *m_helper, flagsPressure, m_params->getPressureAcc(), 0, 0, 0, 1e-4, 1.5, true, 3, false, false, true);
		}
		else {
			if (bivariate) copyVelIntoSrc(velUpdate, flags, m_params->getHeightOfSrc());
			solvePressure(velUpdate, *m_helper, flags, m_params->getPressureAcc(), 0, 0, 0, 1e-4, 1.5, true, 3, false, false, true);
		}  
		if (velPredict) velUpdate.sub(*velPredict);
	}

	KERNEL(idx) void setDenU(const Grid<Real>& denTarget, const Grid<Real>& denPredicut, VectorX& z, const FlagGrid& vh, const int NVel) {
		if (vh(idx) >= 0) {
			z[NVel + vh(idx)] = denTarget(idx)- denPredicut(idx);
		}
	}

	int DivFreeOF::solve(MACGrid& velUpdate, const MACGrid* velPredict, FlagGrid& flags, const Grid<Real>* denPredict, const Grid<Real>* imgs, Tomography* tomo, const Grid<Real>* denTarget) {
		myClock::time_point startTime0 = myClock::now(); 
		velUpdate.clear();
		Grid<Real>* initialDensity;
		if (m_params->isBivariate() && !m_params->shouldUseDenTarget()) initialDensity = new Grid<Real>(velUpdate.getParent());

		// start PD solve
		int iter=0;
		myClock::time_point startTime; 
		Real times[] = { 0.,0.,0. };		

		// adapt flags for pressure solve
		FlagGrid flagsPressure(flags);
		if (m_params->isBivariate()) setSrcToObs(flagsPressure, 0);

		for (iter = 0; iter < m_params->getMaxIter(); iter++) {
			//for (iter = 0; iter < 100; iter++) {
			startTime = myClock::now();
			// *** x-update ***
			// takes its time (res=60, 250ms)
			// first use x+sigma*y as guess, then after enough convergence, use z
			if (iter<7) m_x += m_params->getSigma() * (m_y - m_cg.solveWithGuess(m_x + m_params->getSigma()*m_y - m_b, m_x + m_params->getSigma()*m_y));
			else m_x += m_params->getSigma() * (m_y - m_cg.solveWithGuess(m_x + m_params->getSigma()*m_y - m_b, m_z));
			times[0] += diffClock(myClock::now(), startTime);

			debMsg(iter << ": done with x-update after " << m_cg.iterations() << " CG iterations.", 1);

			// *** z-update ***
			m_z -= m_params->getTau()*m_x;

			// den part if bivariate solve; takes its time (res=60, 110ms)
			if (m_params->isBivariate()) {
				if (m_params->shouldUseDenTarget()) {
					setDenU(*denTarget, *denPredict, m_z, *m_vhDen, m_params->getNumOfVelCmp());
				}
				else {
					// copy values of m_z to initialDensity
					DivFreeOFNS::setInitialDensityClearDenU(*initialDensity, *denPredict, *m_helper, *m_vhDen, m_params->getNumOfVelCmp(), *tomo, m_z);
					tomo->solve(*m_helper, *imgs, initialDensity, &flags);
					DivFreeOFNS::updateDenUpdate(*m_helper, m_z, *denPredict, *m_vhDen, m_params->getNumOfVelCmp(), *tomo);
				}
				times[2] += diffClock(myClock::now(), startTime);
			}

			// vel part
			DivFreeOFNS::copyVecToGrid(velUpdate, *m_vh, m_z, *m_params);
			makeDivFree(velUpdate, velPredict, flags, flagsPressure, denPredict, m_params->getInflowValue(), true, m_params->isBivariate());
			DivFreeOFNS::copyGridToVec(velUpdate, *m_vh, m_z, *m_params);
			times[1] += diffClock(myClock::now(), startTime);

			// *** adapt acc and tolerance ***
			if (iter == (int)ceil(0.33* m_params->getMaxIter()) || iter == (int)ceil(0.66* m_params->getMaxIter()) || iter == m_params->getMaxIter() - 2) {
				m_params->scalePressureAcc(0.2);
				if (m_params->getPressureAcc() < 1e-3) m_params->setPressureAcc(1e-3);
			}

			// *** y-update ***
			m_y = m_z + m_params->getTheta()*(m_z - m_zPrev);
			m_zPrev = m_z;
		}
		
		makeDivFree(velUpdate, velPredict, flags, flagsPressure, denPredict, m_params->getInflowValue(), false, m_params->isBivariate());

		if (m_params->isBivariate()) {
			TomographyNS::copyVecToGridZero(*m_helper, *m_vhDen, m_z, m_params->getNumOfVelCmp());
			// char buff[1000]; snprintf(buff, sizeof(buff), (m_params->getPath() + "densityUpdTomo_%03d_%06d.ppm").c_str(), flags.getSizeX(), m_params->getT());
			// projectPpmFull(*m_helper, buff);
		}

		Real maxDen = (m_params->isBivariate()) ? m_helper->getMaxAbs() : 0;
		printf("t=%03d, velMx=%.3f, denMx=%.3f, solve() took %07.f ms, x-upd: %07.f ms, z-upd vel: %07.f ms, z-upd den: %07.f ms \n", m_params->getT(), velUpdate.getMaxAbs(), maxDen, diffClock(myClock::now(), startTime0), times[0], times[1], times[2]);

		if (m_params->isBivariate() && !m_params->shouldUseDenTarget()) delete initialDensity;
		return iter;
	}

}