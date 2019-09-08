/******************************************************************************
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * @author: Marie-Lena Eckert http://marielenaeckert.com/
 * 
 * Ray for casting rays from pixels through volume
 *
 ******************************************************************************/
#ifndef RAY
#define RAY

#include <omp.h>
#include "manta.h"
#include "commonkernels.h"
#include <Eigen/Sparse>

namespace Manta {

	typedef Vector3D<double> Vec3d;

	inline int getWeights(Real(&weights)[8], const Vec3d& currentPos, const Vec3i& s, const int strideZ) {
		Vec3d fractionsLow, fractionsUp;
		// shift current pos to corner of cell where index starts
		Vec3d posCorner = currentPos - Vec3d(0.5);
		Vec3i idxCorner = Vec3i((int)posCorner.x, (int)posCorner.y, (int)posCorner.z);

		// could parallelize but results in more overhead than speedup
		for (int c = 0; c < 3; c++) {
			if (posCorner[c] < 0.) { idxCorner[c] = 0; fractionsLow[c] = 1.; fractionsUp[c] = 0.; }
			else if (idxCorner[c] >= s[c] - 1) { idxCorner[c] = s[c] - 2; fractionsLow[c] = 0.; fractionsUp[c] = 1.; }
			else { fractionsUp[c] = posCorner[c] - (double)idxCorner[c]; fractionsLow[c] = 1. - fractionsUp[c]; }
		}

		weights[0] = (Real)(fractionsLow.y * fractionsLow.x * fractionsLow.z);
		weights[1] = (Real)(fractionsUp.y * fractionsLow.x * fractionsLow.z);
		weights[2] = (Real)(fractionsLow.y * fractionsUp.x * fractionsLow.z);
		weights[3] = (Real)(fractionsUp.y * fractionsUp.x * fractionsLow.z);
		weights[4] = (Real)(fractionsLow.y * fractionsUp.z * fractionsLow.x);
		weights[5] = (Real)(fractionsUp.y * fractionsUp.z * fractionsLow.x);
		weights[6] = (Real)(fractionsLow.y * fractionsUp.z * fractionsUp.x);
		weights[7] = (Real)(fractionsUp.y * fractionsUp.z * fractionsUp.x);

		return idxCorner.x + s.x * idxCorner.y + strideZ * idxCorner.z;
	}


	class Ray {

	private:
		static const double factorY;       // height = factorY*width
		static const double thresh;	       // thresh for comparing two Vec3 of doubles
		static const double markerWidth;   // double world width of marker, in meters
		static const Vec3d volOffset;      // offset of volume (where volume starts in world coo)
		static const Vec3d volSize;        // size of volume in meters

		// position and direction of ray
		Vec3d m_dir;
		Vec3d m_volEntry;
		Vec3d m_volExit;

		// pixel coordinates
		unsigned short m_x;
		unsigned short m_y;
		unsigned short m_a;

		// check each component of vectors if they are almost equal (thresh)
		inline bool pointsEqual(const Vec3d& a, const Vec3d& b) const;

		// check if point is still in volume cube
		inline bool pointIsInBounds(const Vec3d& p) const;

		// calculate intersection point between line and plane, return Vec3(-1.) if not inside cube (=volume)
		Vec3d iLinePlane(const Vec3d& pos, const Vec3d& p0, const Vec3d& n) const;

		// check intersection points of ray with volume and determine entry (and exit point)
		bool findEntry(const Vec3d& pos);

		// get max number of steps per ray
		inline int getMaxIter(const Vec3d& distance, const Vec3d& stepDir) const;

	public:

		static const double volWidthMeter; // double world size of reconstructed volume

		Ray(const unsigned short x, const unsigned short y, const unsigned short a);

		Ray(const Vec3d& start, const Vec3d& dir, const bool switchXY, const unsigned short x, const unsigned short y, const unsigned short a);

		bool init(const Vec3d& start, const Vec3d& dir, const bool switchXY);

		// get index of pixel, also considering multiple angles
		int getIndex(const unsigned short width, const unsigned short height) const;

		// add up voxel values while tracing, write to separate variable for each ray (pixel value)
		Real traceSum(const Grid<Real>& density, const double stepSize, const FlagGrid* flags) const;

		// setup visual hull for volume while tracing, read and write same variables (entries in vh), but order does not matter
		void traceVH(FlagGrid& vh, const bool blackPixel, const double stepSize) const;

		// fill matrix entries through triplets vector while tracing, write to same variable (triplets) but order does not matter!
		void traceEntries(std::vector<Eigen::Triplet<Real>>& triplets, const FlagGrid& vh, const int idxPixel, const double stepSize, const Real weightAngle) const;

	}; // class

} // namespace

#endif
