/******************************************************************************
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * @author: Marie-Lena Eckert http://marielenaeckert.com/
 *
 * Image with rays casted through volume and concrete pixel values
 *
 ******************************************************************************/

#ifndef IMAGE
#define IMAGE

#include <vector>
#include <iostream>
#include <fstream>
#include "simpleimage.h"
#include "ray.h"
#include "pdSolve.h"

using namespace Eigen;

namespace Manta {

	// -------------------- some helper functions --------------------

	namespace ImageNS {

		enum relation { inside, below, nextToOrBelow };

		inline bool inRelationToShape(const Vec3& pos, const ShapeDetails& shapeLimit, const relation& rel);

		// save images as ppm
		void saveAsPpm(const Grid<Real>& imgsGrid, const std::string& filename);

		void visualizeVH(const FlagGrid& vh, const std::string& filename, const Real scale);

		void saveVHAsGrid(const FlagGrid& vh, const std::string& filename);
	}


	// -------------------- image class --------------------

	PYTHON()
	class Image : public PbClass {

	private:
		unsigned short m_width;
		unsigned short m_height;
		unsigned short m_camCount;
		std::vector<Ray> m_rays;
		std::vector<int> m_raysPerImg;

		// -------------------- some helper methods --------------------

		inline void assertImgSize(const Grid<Real>& imgs) const;
		std::string getBinaryFileName(const std::string& filename, const unsigned short n) const;
		std::string getTxtFilename(const std::string& filename, const unsigned short n) const;
		
	public:
		// filename is required to look like this: ".../mlVid13/calib/%i_rays.txt"
		PYTHON() Image(FluidSolver *parent, const int width, const int height, const int camCount, const std::string& filename, const bool switchXY, const bool orthographic = false);

		// render images from density volume
		PYTHON() void render(Grid<Real>& imgs, const Grid<Real>& density, const double stepSize, const std::string& filenameUni = "", const std::string& filenamePpm = "", const bool sub = false, const FlagGrid* flags = 0) const;

		// setup visual hull for 3D volume based on images and camera rays
		void setupVH(FlagGrid& vh, FlagGrid& vhP, const FlagGrid& flags, const Grid<Real>& imgs, TomoParams& params, const std::string& filename, const Grid<Real>* densityMask) const;
		
		// setup tomography matrix
		void setupMatrix(FlagGrid& vh, FlagGrid& vhP, SparseMatrix<Real, RowMajor>& P, SparseMatrix<Real, RowMajor>& P_single, SparseMatrix<Real, ColMajor>& m_reg, const FlagGrid& flags, const Grid<Real>& imgs, TomoParams& params, const bool usecg, const Grid<Real>* densityMask) const;
		
		unsigned short getCamCount() const;

	}; // class

} // namespace

#endif