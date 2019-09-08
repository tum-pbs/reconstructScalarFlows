/******************************************************************************
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * @author: Marie-Lena Eckert http://marielenaeckert.com/
 *
 * data structures and functions useful to both tomography and divFreeOF PD solves
 *
******************************************************************************/

#include "pdSolve.h"

Eigen::SparseMatrix<Real, Eigen::ColMajor>* m_regLSCG = nullptr;
Eigen::Matrix<Real, -1, 1>* m_rhsUpdateLSCG = nullptr;

using namespace std;

namespace Manta {

	/* ShapeDetails */
	ShapeDetails::ShapeDetails(FluidSolver *parent, const int shape, const Vec3 center, const Vec3 vec, const Real radius) : PbClass(parent), shape(shape), center(center), vec(vec), radius(radius) {}

	int ShapeDetails::getHeightOfSrc() const {
		return int(floor(center.y + vec.y));
	}

	/* PDParams */
	PDParams::PDParams(FluidSolver *parent, const Real sigma, const Real tau, const Real theta, const int maxIter) : PbClass(parent), sigma(sigma), tau(tau), theta(theta), maxIter(maxIter) {}
	
	/* RegWeightsTomo */
	RegWeightsTomo::RegWeightsTomo(FluidSolver *parent, const Real smooth, const Real kinetic, const Real smoothInflow, const Real kineticInflow) : PbClass(parent), smooth(smooth), kinetic(kinetic), smoothInflow(smoothInflow), kineticInflow(kineticInflow) {}

	void RegWeightsTomo::scaleSmooth(Real factor) {
		smooth *= factor;
	}

	void RegWeightsTomo::scaleKinetic(Real factor) {
		kinetic *= factor;
	}

	Real RegWeightsTomo::getSmooth() const {
		return smooth;
	}

	Real RegWeightsTomo::getKinetic() const {
		return kinetic;
	}

	Real RegWeightsTomo::getSmoothInflow() const {
		return smoothInflow;
	}

	Real RegWeightsTomo::getKineticInflow() const {
		return kineticInflow;
	}

	
	/* TomoParams */
	TomoParams::TomoParams(FluidSolver *parent, const int t, const string path, const Real threshVH, const Real threshMask, const Real stepSize, const int minNumCams, const int numPixels, const int numVoxels, const Real angleWeight, const PDParams& pdParams, const RegWeightsTomo& regWeights, const ShapeDetails* shapeLimit)
				: PbClass(parent), t(t), path(path), threshVH(threshVH), threshMask(threshMask), stepSize(stepSize), minNumCams(minNumCams), numPixels(numPixels), numVoxels(numVoxels), angleWeight(angleWeight), pdParams(pdParams), regWeights(regWeights), shapeLimit(shapeLimit) {}
	void TomoParams::printParams() const{
		debMsg(t << ", " << path << ", " << threshVH << ", " << threshMask << ", " << stepSize << ", " << numPixels << ", " << numVoxels << ", " << angleWeight
			<< ", " << pdParams.sigma << ", " << pdParams.tau << ", " << pdParams.theta << ", " << pdParams.maxIter << ", " << getSmooth() << ", " << getKinetic() << ", " << getSmoothInflow() << ", " << getKineticInflow() << ", " << shapeLimit, 1);
	}

	void TomoParams::setT(int newT) {
		t = newT;
	}

	void TomoParams::setNumPixels(int newNumPixels) {
		numPixels = newNumPixels;
	}

	void TomoParams::setNumVoxels(int newNumVoxels) {
		numVoxels = newNumVoxels;
	}
	
	void TomoParams::scaleSmooth(Real factor) {
		regWeights.scaleSmooth(factor);
	}

	void TomoParams::scaleKinetic(Real factor) {
		regWeights.scaleKinetic(factor);
	}

	Real TomoParams::getSmooth() const {
		return regWeights.getSmooth();
	}

	Real TomoParams::getKinetic() const {
		return regWeights.getKinetic();
	}

	Real TomoParams::getSmoothInflow() const {
		return regWeights.getSmoothInflow();
	}

	Real TomoParams::getKineticInflow() const {
		return regWeights.getKineticInflow();
	}

	Real TomoParams::getThreshVH() const {
		return threshVH;
	}

	Real TomoParams::getThreshMask() const {
		return threshMask;
	}

	int TomoParams::getT() const {
		return t;
	}

	string TomoParams::getPath() const {
		return path;
	}

	bool TomoParams::shapeLimitspecified() const {
		return shapeLimit;
	}

	int TomoParams::getNumPixels() const {
		return numPixels;
	}

	int TomoParams::getNumVoxels() const {
		return numVoxels;
	}

	bool TomoParams::zeroVoxels() const {
		return numVoxels == 0;
	}

	Real TomoParams::getStepSize() const {
		return stepSize;
	}

	Real TomoParams::getAngleWeight() const {
		return angleWeight;
	}

	Real TomoParams::getSigma() const {
		return pdParams.sigma;
	}

	Real TomoParams::getTau() const {
		return pdParams.tau;
	}

	Real TomoParams::getTheta() const {
		return pdParams.theta;
	}

	int TomoParams::getMaxIter() const {
		return pdParams.maxIter;
	}

	int TomoParams::getMinNumCams() const {
		return minNumCams;
	}

	const ShapeDetails* TomoParams::getShapeLimit() const {
		return shapeLimit;
	}

	/* RegWeightsOF */
	RegWeightsOF::RegWeightsOF(FluidSolver *parent, const Real smooth, const Real kinetic, const Real kineticZ, const bool adaptiveKineticZ, const Real sumZero) : PbClass(parent), smooth(smooth), kinetic(kinetic), sumZero(sumZero) {
		if (kineticZ < 0) {
			RegWeightsOF::kineticZ = kinetic;
			RegWeightsOF::adaptiveKineticZ = false;
		}
		else {
			RegWeightsOF::kineticZ = kineticZ;
			RegWeightsOF::adaptiveKineticZ = adaptiveKineticZ;
		}
	}

	void RegWeightsOF::scaleSmooth(Real factor) {
		smooth *= factor;
	}

	void RegWeightsOF::scaleKinetic(Real factor) {
		kinetic *= factor;
	}

	Real RegWeightsOF::getSmooth() const {
		return smooth;
	}

	Real RegWeightsOF::getKinetic() const {
		return kinetic;
	}

	Real RegWeightsOF::getKineticZ() const {
		return kineticZ;
	}

	bool RegWeightsOF::useAdaptiveKineticZ() const {
		return adaptiveKineticZ;
	}

	Real RegWeightsOF::getSumZero() const {
		return sumZero;
	}

	/* OFParams */
	OFParams::OFParams(FluidSolver *parent, const int t, const int minSize, const bool bivariate, const string path, const Vec3i s, const Vec3i strides, const int N, const int dim, const int heightOfSrc, const bool useDenTarget, const Real inflowValue, const PDParams& pdParams, const RegWeightsOF& regWeights)
				: PbClass(parent), t(t), minSize(minSize), bivariate(bivariate), path(path), pressureAcc(1e-2), s(s), strides(strides), dim(dim), N(N), heightOfSrc(heightOfSrc), numOfVelCmp(0), useDenTarget(useDenTarget), inflowValue(inflowValue), pdParams(pdParams), regWeights(regWeights) {}

	void OFParams::setT(int newT) {
		t = newT;
	}

	void OFParams::setPressureAcc(Real newPressureAcc) {
		pressureAcc = newPressureAcc;
	}

	void OFParams::scalePressureAcc(Real factor) {
		pressureAcc *= factor;
	}

	void OFParams::setNumOfVelCmp(int newNumOfVelCmp) {
		numOfVelCmp = newNumOfVelCmp;
	}

	void OFParams::setHeightOfSrc(int newHeightOfSrc) {
		heightOfSrc = newHeightOfSrc;
	}

	void OFParams::updateSize(Vec3i newSize) {
		s = newSize;
		strides = Vec3i(1, s.x, s.x*s.y);
		N = s.x*s.y*s.z;
	}

	void OFParams::scaleSmooth(Real factor) {
		regWeights.scaleSmooth(factor);
	}
	
	void OFParams::scaleKinetic(Real factor) {
		regWeights.scaleKinetic(factor);
	}

	Real OFParams::getSmooth() const {
		return regWeights.getSmooth();
	}

	Real OFParams::getKinetic() const {
		return regWeights.getKinetic();
	} 

	Real OFParams::getKineticZ() const {
		return regWeights.getKineticZ();
	}

	Real OFParams::getSumZero() const {
		return regWeights.getSumZero();
	}

	bool OFParams::useAdaptiveKinZ() const {
		return regWeights.useAdaptiveKineticZ();
	}

	int OFParams::getT() const {
		return t;
	}

	Real OFParams::getPressureAcc() const {
		return pressureAcc;
	}

	string OFParams::getPath() const {
		return path;
	}

	int OFParams::getHeightOfSrc() const {
		return heightOfSrc;
	}

	int OFParams::getNumOfVelCmp() const {
		return numOfVelCmp;
	}

	bool OFParams::shouldUseDenTarget() const {
		return useDenTarget;
	}

	Real OFParams::getInflowValue() const {
		return inflowValue;
	}

	int OFParams::getDim() const {
		return dim;
	}

	int OFParams::getMaxIter() const {
		return pdParams.maxIter;
	}

	Vec3i OFParams::getStrides() const {
		return strides;
	}

	int OFParams::getStride(int c) const {
		return strides[c];
	}

	int OFParams::getN() const {
		return N;
	}

	int OFParams::getS(int c) const {
		return s[c];
	}

	int OFParams::getMinSize() const {
		return minSize;
	}

	Real OFParams::getTheta() const {
		return pdParams.theta;
	}

	Real OFParams::getSigma() const {
		return pdParams.sigma;
	}

	Real OFParams::getTau() const {
		return pdParams.tau;
	}


	bool OFParams::isBivariate() const {
		return bivariate;
	}

	void OFParams::printParams() const {
		debMsg(t<<", "<<minSize<<", "<<bivariate<<", "<<path<<", "<<pressureAcc<<", "<<s<<", "<<strides<<", "<<dim<<", "<<N<<", "<<heightOfSrc<<", "<<numOfVelCmp 
			<< ", " << pdParams.sigma << ", " << pdParams.tau << ", " << pdParams.theta << ", " << pdParams.maxIter << ", " << getSmooth() << ", " << getKinetic() << ", " << getKineticZ() << ", " << useAdaptiveKinZ() << ", " << getSumZero(),1 );
	}
}