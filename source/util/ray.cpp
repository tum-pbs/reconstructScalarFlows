/******************************************************************************
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * @author: Marie-Lena Eckert http://marielenaeckert.com/
 *
 * Ray for casting rays from pixels through volume
 *
 ******************************************************************************/

#include "ray.h"

namespace Manta {

	const double Ray::factorY = 1.77;
	const double Ray::thresh = 1e-12;
	const double Ray::markerWidth = 0.4909;
	const double Ray::volWidthMeter = 1.*Ray::markerWidth;
	
	// move volume 
	const Vec3d Ray::volOffset = Vec3d(markerWidth / 6., -markerWidth / 11, -markerWidth / 100);// 1 * (markerWidth - volWidthMeter), -markerWidth / 15., 1 * (markerWidth - volWidthMeter));
	const Vec3d Ray::volSize(volWidthMeter, ceil(factorY*volWidthMeter), volWidthMeter);

	inline bool Ray::pointsEqual(const Vec3d& a, const Vec3d& b) const {
		if (fabs(a.x - b.x) > thresh) return false;
		if (fabs(a.y - b.y) > thresh) return false;
		if (fabs(a.z - b.z) > thresh) return false;
		return true;
	}

	// check if point is still in volume cube
	inline bool Ray::pointIsInBounds(const Vec3d& p) const {
		if (p.x < -thresh + volOffset.x || p.x > volOffset.x + volSize.x + thresh
		 || p.y < -thresh + volOffset.y || p.y > volOffset.y + volSize.y + thresh
		 || p.z < -thresh + volOffset.z || p.z > volOffset.z + volSize.z + thresh) {
			return false;
		}
		return true;
	}

	// calculate intersection point between line and plane, return Vec3(-1.) if not inside cube (=volume)
	Vec3d Ray::iLinePlane(const Vec3d& pos, const Vec3d& p0, const Vec3d& n) const {
		// check if ray and plane are parallel
		double divideBy = dot(m_dir, n);
		if (fabs(divideBy) < thresh) return Vec3d(-1.);

		Vec3d p = pos + (dot((p0 - pos), n) / divideBy) * m_dir;

		if (!pointIsInBounds(p)) return Vec3d(-1.);

		return p;
	}

	// check intersection points of ray with volume and determine entry and exit point
	bool Ray::findEntry(const Vec3d& pos) {
		// calculate intersection point for all 6 bounding planes
		Vec3d p[] = {
			iLinePlane(pos, volOffset,			 Vec3d(1.,0.,0.)),
			iLinePlane(pos, volOffset + volSize, Vec3d(1.,0.,0.)),
			iLinePlane(pos, volOffset,			 Vec3d(0.,1.,0.)),
			iLinePlane(pos, volOffset + volSize, Vec3d(0.,1.,0.)),
			iLinePlane(pos, volOffset,			 Vec3d(0.,0.,1.)),
			iLinePlane(pos, volOffset + volSize, Vec3d(0.,0.,1.)) };

		// determine which intersection points are entry and exit point
		int nIntersec = 0;
		double entryDist = std::numeric_limits<double>::max();
		for (int i = 0; i < 6; i++) {
			if (!pointsEqual(p[i], Vec3d(-1.))) {
				nIntersec++;
				double dist = norm(p[i] - pos);
				if (dist < entryDist) {
					entryDist = dist;
					m_volEntry = p[i];
					if (pointsEqual(m_volExit, Vec3d(-1.))) m_volExit = m_volEntry;
				}
				else if (dist - thresh > entryDist) m_volExit = p[i];
			}
		}

		// ensure valid ray-volume intersection

		// standard case: 2 intersections with volume, ensure entry and exit points are valid; if entry == exit, only one intersection took place
		if (nIntersec == 2 && !pointsEqual(m_volEntry, m_volExit)) {
			if (pointsEqual(m_volEntry, Vec3d(-1.)) || pointsEqual(m_volExit, Vec3d(-1.))) {
				errMsg("findEntry: entry and exit point are not valid: " << m_volEntry << ", " << m_volExit);
				return false;
			}
			return true;
		}// check if more than 2 intersections with volume sides and if they're not the same
		else {
			for (int i = 0; i < 6; i++) {
				if (!pointsEqual(p[i], Vec3d(-1.)) && !pointsEqual(p[i], m_volEntry) && !pointsEqual(p[i], m_volExit)) {
					debMsg(p[0] << ", " << p[1] << ", " << p[2] << ", " << p[3] << ", " << p[4] << ", " << p[5], 1);
					errMsg("findEntry: more than two intersections between ray and volume (cube) and points do not match: " << nIntersec << ", " << m_volEntry << ", " << m_volExit << ", " << p[i]);
				}
			}
			return false;
		}
	}

	inline int Ray::getMaxIter(const Vec3d& distance, const Vec3d& stepDir) const {
		return floor(distance.x / stepDir.x);
	}

	Ray::Ray(const unsigned short x, const unsigned short y, const unsigned short a) : m_dir(Vec3d(-1)), m_volEntry(Vec3d(-1)), m_volExit(Vec3d(-1)), m_x(x), m_y(y), m_a(a) {}

	Ray::Ray(const Vec3d& start, const Vec3d& dir, const bool switchXY, const unsigned short x, const unsigned short y, const unsigned short a) : m_volEntry(Vec3d(-1)), m_volExit(Vec3d(-1)), m_x(x), m_y(y), m_a(a) {
		init(start, dir, switchXY);
	}

	bool Ray::init(const Vec3d& start, const Vec3d& dir, const bool switchXY) {
		if (switchXY) m_dir = Vec3d(dir.y, dir.x, dir.z);
		else m_dir = dir;

		// pos in world space, in meters
		Vec3d pos = - 50. * m_dir; // move pos away from front marker plane to ensure pos is outside the reconstructed volume; could also use 100*dir
		if (switchXY)  pos += Vec3d(start.y, start.x, start.z);
		else pos += start;

		// find entry and exit point of ray through volume
		// does not depend on grid size, takes place in world coordinates
		if (findEntry(pos)) {
			m_volEntry -= volOffset;
			m_volExit -= volOffset;
			return true;
		}
		else return false;
	}

	int Ray::getIndex(const unsigned short width, const unsigned short height) const {
		return m_x + m_y*width + m_a*width*height;
	}

	// extensively used methods, should be as fast as possible

	// add up voxel values while tracing, used once in each inner PD loop
	Real Ray::traceSum(const Grid<Real>& density, const double stepSize, const FlagGrid* flags) const {
		const Vec3i s = density.getSize();
		const int strideZ = s.x*s.y;
		const int neighbours[8] = { 0 , s.x, 1, 1 + s.x, strideZ, s.x + strideZ, 1 + strideZ, 1 + s.x + strideZ };

		const Vec3d stepDir = stepSize*m_dir;
		const double realToVol = s.x / volWidthMeter;
		const Vec3d startPos = m_volEntry*realToVol;
		const int maxIter = getMaxIter((m_volExit - m_volEntry)*realToVol, stepDir);

		Real weights[8];
		Real pixelVal = 0;
		for (int i = 0; i < maxIter; i++) {
			Vec3d currentPos = startPos + i * stepDir;
			int idx = getWeights(weights, currentPos, s, strideZ);
			for (int ind = 0; ind < 8; ind++) {
				// don't care if weight is almost zero, since this check is more expensive than the multiplication with 0
				if (!flags || !flags->isSrc(idx + neighbours[ind])) pixelVal += weights[ind] * density(idx + neighbours[ind]);
			}
		}

		// final scaling
		return (stepSize / s.x) * pixelVal;
	}

	// setup visual hull for volume while tracing, only used once per time step
	void Ray::traceVH(FlagGrid& vh, const bool blackPixel, const double stepSize) const {
		const Vec3i s = vh.getSize();
		const Real mult = (stepSize / s.x); 
		const int strideZ = s.x*s.y;
		const int neighbours[8] = { 0 , s.x, 1, 1 + s.x, strideZ, s.x + strideZ, 1 + strideZ, 1 + s.x + strideZ };

		const Vec3d stepDir = stepSize*m_dir;
		const double realToVol = s.x / volWidthMeter;
		const Vec3d startPos = m_volEntry*realToVol;
		const int maxIter = getMaxIter((m_volExit-m_volEntry)*realToVol, stepDir);

		Real weights[8];

		// could parallelize more, but only used once and shared vh (read and write)
		for (int i = 0; i < maxIter; i++) {
			Vec3d currentPos = startPos + i * stepDir;
			int idx = getWeights(weights, currentPos, s, strideZ);
			for (int ind = 0; ind < 8; ind++) {
				// read and write vh
				// if already -1 or weight almost 0, don't process
				if (fabs(mult * weights[ind]) < 1e-4 || vh(idx + neighbours[ind]) == -1) continue;
				if (blackPixel) {
					// don't remove voxel from visual hull if weight very small, could have some value but still seen as 'black'
					if (fabs(weights[ind]) > 1e-1) vh(idx + neighbours[ind]) = -1;
				}
				else vh(idx + neighbours[ind])++;
			}
		}
	}

	// fill matrix entries through triplets vector while tracing, used in each inner PD loop in every CG loop multiple times
	void Ray::traceEntries(std::vector<Eigen::Triplet<Real>>& triplets, const FlagGrid& vh, const int idxPixel, const double stepSize, const Real weightAngle) const {
		const Vec3i s = vh.getSize();
		const Real mult = (stepSize / s.x);
		const int strideZ = s.x*s.y;
		const int neighbours[8] = { 0 , s.x, 1, 1 + s.x, strideZ, s.x + strideZ, 1 + strideZ, 1 + s.x + strideZ };

		const Vec3d stepDir = stepSize*m_dir;
		const double realToVol = s.x / volWidthMeter;
		const Vec3d startPos = m_volEntry*realToVol;
		const int maxIter = getMaxIter((m_volExit - m_volEntry)*realToVol, stepDir);

		triplets.reserve(4 * maxIter); // max size is 8*maxIter - visual hull - small weight

		Real weights[8];

		// loop is not iteration independent! pushes back to same vector
		for (int i = 0; i < maxIter; i++) {
			Vec3d currentPos = startPos + i * stepDir;
			int idx = getWeights(weights, currentPos, s, strideZ);
			for (int ind = 0; ind < 8; ind++) {
				// if weight almost 0, don't process
				if (fabs(mult * weights[ind]) < 1e-4) continue;
				int idxVoxel = vh(idx + neighbours[ind]);
				// if already -1, don't process
				if (idxVoxel == -1) continue;

				if (m_a == 0) triplets.push_back(Eigen::Triplet<Real>(idxPixel, idxVoxel, mult * weights[ind]));
				else triplets.push_back(Eigen::Triplet<Real>(idxPixel, idxVoxel, weightAngle * mult * weights[ind]));
			}
		}
	}

} // namespace
