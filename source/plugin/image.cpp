/******************************************************************************
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * @author: Marie-Lena Eckert http://marielenaeckert.com/
 *
 * Image with rays casted through volume and concrete pixel values
 *
******************************************************************************/

#include "image.h"

namespace Manta {

	// -------------------- Kernel functions --------------------

	namespace ImageNS {
		// copy images from uni grid to ppm images
		KERNEL() void copyUniToPpm(const Grid<Real>& imgsGrid, std::vector<SimpleImage>& imgsPpm) {
			imgsPpm[k](i, j) = imgsGrid(i, j, k);
		}

		// save flag grid as real grid in order to visualize it nicely
		KERNEL(idx) void flagToReal(Grid<Real>& density, const FlagGrid& vh) {
			if (vh(idx) == -1) density(idx) = 0;
			else density(idx) = 1;
		}

		// returns false if pos is not inside/below/nextToOrBelow shape specified in shapeLimit
		inline bool inRelationToShape(const Vec3& pos, const ShapeDetails& shapeLimit, const relation& rel) {
			if (shapeLimit.shape == 0) { // cylinder 
				Vec3 axis = shapeLimit.vec;
				Real length = normalize(axis); // axis is now normalized
				Real z = dot(pos - shapeLimit.center, axis);

				switch (rel)
				{
				case inside:
					if (fabs(z) > length) return false;
					else {
						Real r2 = normSquare(pos - shapeLimit.center) - square(z);
						return r2 < square(shapeLimit.radius);
					}
				default:
					errMsg("inRelationToShape: only relation inside is implemented for source shape 0 (=cylinder)");
					return true;
				}
			}
			else if (shapeLimit.shape == 1) { // sphere
				switch (rel)
				{
				case inside:
					return normSquare((pos - shapeLimit.center)) <= shapeLimit.radius * shapeLimit.radius;
				default:
					errMsg("inRelationToShape: only relation inside is implemented for source shape 1 (=sphere)");
					return true;
				}
			}
			else if (shapeLimit.shape == 2) { // box
				Vec3 mP0 = shapeLimit.center - shapeLimit.vec;
				Vec3 mP1 = shapeLimit.center + shapeLimit.vec;
				switch (rel)
				{
				case inside:
					return pos.x >= mP0.x && pos.y >= mP0.y && pos.z >= mP0.z && pos.x <= mP1.x && pos.y <= mP1.y && pos.z <= mP1.z;
				case below:
					return pos.x >= mP0.x && pos.y < mP0.y && pos.z >= mP0.z && pos.x <= mP1.x && pos.z <= mP1.z;
				case nextToOrBelow: // true if either below box or next to box
					return pos.y < mP0.y || (pos.y <= mP1.y && (pos.x < mP0.x || pos.x > mP1.x || pos.z < mP0.z || pos.z > mP1.z));
				default:
					return true;
				}
			}
			else return true;
		}

		// postprocess vh for pixels
		KERNEL(single, idx, reduce = +) returns(int idxPixel = 0)
		int finalizePixelVH(FlagGrid& vhP) {
			if(vhP(idx)==1) vhP(idx) = idxPixel++;
		}

		// create an index for all seen, non-black voxels considering shape, bnds, mask
		KERNEL(single, reduce = +) returns(int idxVoxel = 0)
		int finalizeVoxelVH(FlagGrid& vh, const std::vector<FlagGrid*>& vhs, const FlagGrid& flags, const Grid<Real>* densityMask, const Real threshMask, const ShapeDetails* shapeLimit, const relation& rel, const int minNumCams) {
			int numSeeBlack = 0;
			int numOfCamsSeeing = 0;
			for (auto vhA : vhs) {
				if ((*vhA)(i, j, k) > 0) numOfCamsSeeing++; 
				else if ((*vhA)(i, j, k) == -1) numSeeBlack++;
			}

			if (numSeeBlack>=1 || numOfCamsSeeing < minNumCams || !flags.isFluid(i, j, k) || flags.isSrc(i, j, k) || (threshMask > 0 && densityMask && fabs((*densityMask)(i, j, k)) < threshMask)
				|| (shapeLimit!=nullptr && shapeLimit->shape >= 0 && shapeLimit->shape <= 2 && !inRelationToShape(Vec3(i + 0.5, j + 0.5, k + 0.5), *shapeLimit, rel))) {
				vh(i, j, k) = -1;
			}
			else {
				// now create index 
				vh(i, j, k) = idxVoxel++;
			}
		}

		// -------------------- helper functions --------------------

		// save images as ppm
		PYTHON() void saveAsPpm(const Grid<Real>& imgsGrid, const std::string& filename) {
			// init simple images, one for each camera
			std::vector<SimpleImage> imgsPpm(imgsGrid.getSizeZ());
			for (auto it = imgsPpm.begin(); it != imgsPpm.end(); ++it) (*it).init(imgsGrid.getSizeX(), imgsGrid.getSizeY());

			// copy values
			copyUniToPpm(imgsGrid, imgsPpm);

			// write to disk
			char buff[1000];
			for (int i = 0; i < imgsPpm.size(); i++) {
				snprintf(buff, sizeof(buff), filename.c_str(), i + 1);
				imgsPpm[i].writePpm(buff);
			}
		}

		// visualize visual hull as real grid projected to ppm
		PYTHON() void visualizeVH(const FlagGrid& vh, const std::string& filename, const Real scale) {
			Grid<Real> density(vh.getParent());
			flagToReal(density, vh);
			projectPpmFull(density, filename, 0, scale);
		}

		// visualize visual hull as real grid projected to ppm
		void saveVHAsGrid(const FlagGrid& vh, const std::string& filename) {
			Grid<Real> density(vh.getParent());
			flagToReal(density, vh);
			density.save(filename);
		}
	}


	// -------------------- class method definitions --------------------

	// some helper methods 

	// ensure that class members match grid size
	inline void Image::assertImgSize(const Grid<Real>& imgs) const {
		if (imgs.getSizeX() != m_width || imgs.getSizeY() != m_height || imgs.getSizeZ() != m_camCount) {
			errMsg("Error in render image, images size does not match: " << imgs.getSize() << " vs. (" << m_width << ", " << m_height << ", " << m_camCount << ")");
		}
	}

	std::string Image::getBinaryFileName(const std::string& filename, const unsigned short n) const {
		char buff[1000]; // needed for string formatting for filename
		snprintf(buff, sizeof(buff), (filename.substr(0, filename.length() - 4) + "_%d_%d.bin").c_str(), n, m_width, m_height);
		return buff;
	}

	std::string Image::getTxtFilename(const std::string& filename, const unsigned short n) const {
		char buff[1000]; // needed for string formatting for filename
		snprintf(buff, sizeof(buff), filename.c_str(), n);
		return buff;
	}

	// main methods

	// constructor: initialize all rays, calculate entry point
	Image::Image(FluidSolver *parent, const int width, const int height, const int camCount, const std::string& filename, const bool switchXY, const bool orthographic) : PbClass(parent), m_width(width), m_height(height), m_camCount(camCount)
	{
		int estimCntRays = m_width*m_height*m_camCount;

		// for each camera image, read ray file or setup rays manually if orthographic
		m_raysPerImg = std::vector<int>(m_camCount, 0);
		for (unsigned short a = 0; a < m_camCount; a++) {
			if (orthographic) {
				m_rays.reserve(estimCntRays);
				Real scalePixels = Ray::volWidthMeter / m_width;
				Real shiftPixels = Ray::volWidthMeter / 2.;
				for (unsigned short x = 0; x < m_width; x++) {
					for (unsigned short y = 0; y < m_height; y++) {
						double xPixel = x * scalePixels + shiftPixels;
						double yPixel = y * scalePixels + shiftPixels;
						Ray ray(x,y,a);
						if (a == 0) {
							if (ray.init(Vec3d(xPixel, yPixel, -10), Vec3d(0, 0, 1), switchXY)) {
								m_rays.push_back(ray);
								m_raysPerImg[a]++;
							}
						}
						else if (a == 1) {
							if (ray.init(Vec3d(-10, yPixel, xPixel), Vec3d(1, 0, 0), switchXY)) {
								m_rays.push_back(ray);
								m_raysPerImg[a]++;
							}
						}
						else {
							if (ray.init(Vec3d(xPixel, -10, yPixel), Vec3d(0, 1, 0), switchXY)) {
								m_rays.push_back(ray);
								m_raysPerImg[a]++;
							}
						}
					}
				}
			}
			else {
				m_rays.reserve(estimCntRays*0.6);
				std::string binaryFilename = getBinaryFileName(filename, a + 1);
				std::fstream fileBinary(binaryFilename, std::ios::in | std::ios::binary);
				// if binary file exists, read in
				if (fileBinary.good()) {
					while (!fileBinary.eof()) {
						unsigned short x, y;
						fileBinary.read(reinterpret_cast<char*>(&x), sizeof(unsigned short));
						fileBinary.read(reinterpret_cast<char*>(&y), sizeof(unsigned short));
						Vec3d start, dir;
						for (int c = 0; c < 3; c++) {
							fileBinary.read(reinterpret_cast<char*>(&start[c]), sizeof(double));
							fileBinary.read(reinterpret_cast<char*>(&dir[c]), sizeof(double));
						}
						Ray ray(x,  y , a);
						if (ray.init(start, dir, switchXY)) {
							m_rays.push_back(ray);
							m_raysPerImg[a]++;
						}
					}
				}
				// otherwise write binary file
				else {
					fileBinary.open(binaryFilename, std::ios::out | std::ios::binary);
					std::string txtFilename = getTxtFilename(filename, a + 1);
					std::ifstream rayfile(txtFilename);

					unsigned short sX = 0, sY = 0;
					rayfile >> sX >> sY;
					debMsg("Reading ray file " << txtFilename << " with w=" << sX << " and h=" << sY << " -> scale to w=" << m_width << " h=" << m_height, 1);

					// number of rays mapped to an output pixel
					unsigned int* cntRays = new unsigned int[m_width*m_height](); // initialized to zero
					Vec3d* starts = new Vec3d[m_width*m_height];                      // not initialized
					Vec3d* dirs = new Vec3d[m_width*m_height];                        // not initialized

					// now read the data
					unsigned short x, y;
					double startX, startY, dirX, dirY, dirZ;
					while (rayfile >> x >> y >> startX >> startY >> dirX >> dirY >> dirZ) {
						// calculate the point in scaled space
						x = (x * width) / sX;
						y = (y * height) / sY;

						// check if pixel is out of bounds
						if (x < 0 || x >= m_width || y < 0 || y >= m_height) {
							std::cerr << "Ignoring out-of-bounds point with x=" << x << " y=" << y << std::endl;
							continue;
						}

						// index of pixel
						unsigned int idx = x + m_width*y;
						if (cntRays[idx] == 0) {
							// set start position and direction vector
							starts[idx] = Vec3d(startX, startY, 0.);
							dirs[idx] = Vec3d(dirX, dirY, dirZ);
						}
						else {
							// sum up different start positions and direction vectors to build the average later
							starts[idx] += Vec3d(startX, startY, 0.);
							dirs[idx] += Vec3d(dirX, dirY, dirZ);
						}
						cntRays[idx]++;
					}

					for (x = 0; x < m_width; x++) {
						for (y = 0; y < m_height; y++) {
							// index of pixel
							unsigned int idx = x + m_width*y;
							if (cntRays[idx] == 0) continue;

							Vec3d start = starts[idx] / cntRays[idx];
							Vec3d dir = dirs[idx] / norm(dirs[idx]);

							fileBinary.write(reinterpret_cast<char*>(&x), sizeof(unsigned short));
							fileBinary.write(reinterpret_cast<char*>(&y), sizeof(unsigned short));
							for (int c = 0; c < 3; c++) {
								fileBinary.write(reinterpret_cast<char*>(&start[c]), sizeof(double));
								fileBinary.write(reinterpret_cast<char*>(&dir[c]), sizeof(double));
							}

							Ray ray(x, y, a);
							if (ray.init(start, dir, switchXY)) {
								m_rays.push_back(ray);
								m_raysPerImg[a]++;
							}
						}
					}

					delete[] cntRays;
					delete[] starts;
					delete[] dirs;
				}
				fileBinary.close();
			}
		}

		m_rays.shrink_to_fit();
		for (int a = 1; a < m_camCount; a++) m_raysPerImg[a] += m_raysPerImg[a - 1];
	}

	// render images from density volume
	void Image::render(Grid<Real>& imgs, const Grid<Real>& density, const double stepSize, const std::string& filenameUni, const std::string& filenamePpm, const bool sub, const FlagGrid* flags) const {
		assertImgSize(imgs);
		if (!sub) imgs.clear();

#pragma omp parallel for
		for (int i = 0; i<m_rays.size(); i++) {
			if(sub) imgs(m_rays[i].getIndex(m_width, m_height)) -= m_rays[i].traceSum(density, stepSize, flags);
			else imgs(m_rays[i].getIndex(m_width, m_height)) = m_rays[i].traceSum(density, stepSize, flags);
		}

		if (filenameUni != "") imgs.save(filenameUni);
		if (filenamePpm != "") ImageNS::saveAsPpm(imgs, filenamePpm);
	}

	inline bool isNeighborhoodBlack(const Grid<Real>& imgs, const int idxCurrentPixel, const int width, const int height, const int wh, const Real threshVH, const int resX) {
		int x = idxCurrentPixel % width;
		int y = ((idxCurrentPixel - x) % wh) / width;
		int z = idxCurrentPixel / wh;
		int kernelSize = 31;//width/resX*5 + 1;// 13;
		if (kernelSize % 2 == 0) kernelSize++;
		int center = kernelSize / 2;
		for (int j = 0; j < kernelSize; j++) {
			for (int i = 0; i < kernelSize; i++) {
				int xx = x + i - center;
				int yy = y + j - center;
				if (xx >= 0 && yy >= 0 && xx <= width - 1 && yy <= height - 1) {
					if (fabs(imgs(xx, yy, z)) > threshVH) return false;
				}
			}
		}
		return true;
	}

	inline void processRayVH(std::vector<FlagGrid*> vhs, FlagGrid& vhP, const Grid<Real>& imgs, const Ray& ray, const int width, const int height, const int wh, const Real& threshVH, const Real& stepSize) {
		int idxCurrentPixel = ray.getIndex(width, height);
		// check neighbours as well in order to avoid resolution issues (e.g. 4 pixels for same voxels, some are zero, some aren't)
		bool isBlackPixel = (fabs(imgs(idxCurrentPixel)) < threshVH) && isNeighborhoodBlack(imgs, idxCurrentPixel, width, height, wh, threshVH, vhs[0]->getSizeX());

		// create new index for pixels
		if (isBlackPixel) vhP(idxCurrentPixel) = -1;
		else vhP(idxCurrentPixel) = 1;

		// trace through volume
		ray.traceVH(*vhs[idxCurrentPixel / wh], isBlackPixel, stepSize);
	}

	// setup visual hull for 3D volume based on images and camera rays
	void Image::setupVH(FlagGrid& vh, FlagGrid& vhP, const FlagGrid& flags, const Grid<Real>& imgs, TomoParams& params, const std::string& filename, const Grid<Real>* densityMask) const {
		assertImgSize(imgs);
		
		std::vector<FlagGrid*> vhs; vhs.reserve(imgs.getSizeZ());
		for (int a = 0; a < imgs.getSizeZ(); a++) {
			vhs.push_back(new FlagGrid(vh.getParent()));
			vhs[a]->setConst(0);
		}
		if (vhs.size() == 0) errMsg("Error in setupVH: zero camera angles.");

		int idxCurrentPixel = 0;
		int wh = m_height*m_width;
#pragma omp parallel for
		for (int a = 0; a < m_camCount; a++) {
			int start = (a == 0) ? 0 : m_raysPerImg[a - 1];
			for (int r = start; r < m_raysPerImg[a]; r++) {
				processRayVH(vhs, vhP, imgs, m_rays[r], m_width, m_height, wh, params.getThreshVH(), params.getStepSize());
			}
		}
		
		// vhP now contains values of -1 for black pixels and index values for non-black pixels
		// vh now contains values of -1 voxels seen as black, 0 for unseen voxels, >0 for non-black voxels

		// now create an index for all seen, non-black voxels
		// also consider shapeLimit, bnds, mask
		params.setNumVoxels(ImageNS::finalizeVoxelVH(vh, vhs, flags, densityMask, params.getThreshMask(), params.getShapeLimit(), ImageNS::relation::inside, params.getMinNumCams()));
		params.setNumPixels(ImageNS::finalizePixelVH(vhP));

		if(filename!="") ImageNS::visualizeVH(vh, filename, 1);

		for (int a = 0; a < imgs.getSizeZ(); a++) delete vhs[a];
	}

	inline void addRhoRegs(std::vector<Triplet<Real>>& triplets, const FlagGrid& vh, const FlagGrid& flags, const TomoParams& params, const Vec3i& s, const int dim) {
		// add rho and regs
		Vec3i strides = Vec3i(1, s.x, s.x*s.y);
		Real mult;
		for (int k = 0; k < s.z; k++) {
			for (int j = 0; j < s.y; j++) {
				for (int i = 0; i < s.x; i++) {
					int ijk = vh.index(i, j, k);
					int idx = vh(ijk);
					if (idx == -1 || !vh.isInBounds(Vec3i(i,j,k), 1)) continue;

					// smoothness and kinetic energy regularizer if neighbour in visual hull
					mult = 0;
					for (int c = 0; c < dim; c++) {
						int idxLow = vh(ijk - strides[c]);
						int idxUp = vh(ijk + strides[c]);
						if (idxLow != -1) {
							triplets.push_back(Triplet<Real>(idx, idxLow, -params.getSmooth()));
							mult++;
						}
						if (idxUp != -1) {
							triplets.push_back(Triplet<Real>(idx, idxUp, -params.getSmooth()));
							mult++;
						}
					}
					triplets.push_back(Triplet<Real>(idx, idx, params.getSigma() + mult*params.getSmooth() + params.getKinetic()));
				}
			}
		}
	}

	// setup tomography matrix: P, PP+rhoI
	void Image::setupMatrix(FlagGrid& vh, FlagGrid& vhP, SparseMatrix<Real, RowMajor>& P, SparseMatrix<Real, RowMajor>& P_single, SparseMatrix<Real, ColMajor>& m_reg, const FlagGrid& flags, const Grid<Real>& imgs, TomoParams& params, const bool usecg, const Grid<Real>* densityMask) const {
		myClock::time_point startTime0 = myClock::now();
		myClock::time_point startTime = startTime0;

		// create vh in order to know which pixel-voxel entries to save in matrix
		setupVH(vh, vhP, flags, imgs, params, "", densityMask);

		// setup triplet vector
		Vec3i s = vh.getSize();
		// setup regularizing matrix for LSCG solve
		std::vector<Triplet<Real>> triplets;
		triplets.reserve(5 * params.getNumVoxels());
		addRhoRegs(triplets, vh, flags, params, s, 3);
		m_reg.resize(params.getNumVoxels(), params.getNumVoxels());
		m_reg.setFromTriplets(triplets.begin(), triplets.end());

		triplets.clear();
		triplets.shrink_to_fit();
		triplets.reserve((17.)*params.getNumVoxels()*s.z/params.getStepSize());

		// go through rays from non-black pixels and save pixel-voxel entries
		// not parallel since triplets are written at the same time; 
		for (int i = 0; i<m_rays.size(); i++) {
			int idxPixel = vhP(m_rays[i].getIndex(m_width, m_height));
			if (idxPixel == -1) continue;
			m_rays[i].traceEntries(triplets, vh, idxPixel, params.getStepSize(), params.getAngleWeight());
		}

		// resize P to nonzero-pixels and vh-voxels
		P.resize(params.getNumPixels(), params.getNumVoxels());
		P.setFromTriplets(triplets.begin(), triplets.end());

		if (usecg) {
			P_single = P;
			P = P_single.transpose()*P_single + m_reg;
		}
	}

	unsigned short Image::getCamCount() const {
		return m_camCount;
	}

} // namespace
