/*
By downloading, copying, installing or using the software you agree to this
license. If you do not agree to this license, do not download, install,
copy or use the software.

                          License Agreement
               For Open Source Computer Vision Library
                       (3-clause BSD License)

Copyright (C) 2013, OpenCV Foundation, all rights reserved.
Third party copyrights are property of their respective owners.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

  * Redistributions of source code must retain the above copyright notice,
    this list of conditions and the following disclaimer.

  * Redistributions in binary form must reproduce the above copyright notice,
    this list of conditions and the following disclaimer in the documentation
    and/or other materials provided with the distribution.

  * Neither the names of the copyright holders nor the names of the contributors
    may be used to endorse or (setq c-basic-offset 4)promote products derived from this software
    without specific prior written permission.

This software is provided by the copyright holders and contributors "as is" and
any express or implied warranties, including, but not limited to, the implied
warranties of merchantability and fitness for a particular purpose are
disclaimed. In no event shall copyright holders or contributors be liable for
any direct, indirect, incidental, special, exemplary, or consequential damages
(including, but not limited to, procurement of substitute goods or services;
loss of use, data, or profits; or business interruption) however caused
and on any theory of liability, whether in contract, strict liability,
or tort (including negligence or otherwise) arising in any way out of
the use of this software, even if advised of the possibility of such damage.
*/


/*
 * This is based on the the OpenCV example detect_board_charuco
 */

#include <opencv2/highgui.hpp>
#include <opencv2/aruco/charuco.hpp>
#include <opencv2/opencv.hpp>
#include <vector>
#include <iostream>
#include <fstream>
#include <tbb/tbb.h>

using namespace std;
using namespace cv;
using namespace tbb;

/*
 * Definition of a single ray.
 */
struct EstimatedRay
{
    /*
     * Start point
     */
    Vec3f start;

    /*
     * Direction (normalized)
     */
    Vec3f dir;

    /*
     * This is an indicator for the quality / precision expected for this ray.
     *
     * A ray is better if the source points:
     *   - are NOT extrapolated!
     *   - are closer to a reference point
     *   - have a higher Z distance
     *
     * Smaller values are better.
     * Range is from 0 to infinity.
     */
    float rating;

    // Operator for sorting by rating
    bool operator<(const EstimatedRay &other) const
    {
        return rating < other.rating;
    }
};


/*
 * A MappedQuad represents a quad with known world and image positons.
 */
struct MappedQuad
{
    /*
     * The four corner points in the image.
     */
    Point2f corners_image[4];

    /*
     * The two edges of the worldspace rect. This is always axis aligned.
     */
    Point2f ws_min, ws_max; 

    /*
     * Indices of the neighbouring quads.
     * Nonexistant quads are set to MAX_SIZE.
     */
    size_t adjacent_quads[8];

    /*
     * Tranformation from image point to world point.
     * P_world = transform * P_image
     */
    Mat transform;
};

namespace {
const char* about = "Pose estimation using a ChArUco board";
const char* keys  =
        "{w        |       | Number of squares in X direction }"
        "{h        |       | Number of squares in Y direction }"
        "{sl       |       | Square side length (in meters) }"
        "{ml       |       | Marker side length (in meters) }"
        "{d        |       | dictionary: DICT_4X4_50=0, DICT_4X4_100=1, DICT_4X4_250=2,"
        "DICT_4X4_1000=3, DICT_5X5_50=4, DICT_5X5_100=5, DICT_5X5_250=6, DICT_5X5_1000=7, "
        "DICT_6X6_50=8, DICT_6X6_100=9, DICT_6X6_250=10, DICT_6X6_1000=11, DICT_7X7_50=12,"
        "DICT_7X7_100=13, DICT_7X7_250=14, DICT_7X7_1000=15, DICT_ARUCO_ORIGINAL = 16}"
        "{rs       |       | Apply refind strategy }"
        "{dp       |       | File of marker detector parameters }"
        "{v        |       | Input files / images}"
        "{z        | 1     | Depth distance between two steps}"
        "{o        |       | Prefix of the output files}"
        "{ed       | 200   | Max. distance in pixels for extraploation }"
        "{ee       | 0.025 | Max. error in extrapolation }"
        "{em       | 3     | Minimum number of quads in range for extrapolation }";
}

/*
 */
static bool readDetectorParameters(string filename, Ptr<aruco::DetectorParameters> &params) {
    FileStorage fs(filename, FileStorage::READ);
    if(!fs.isOpened())
        return false;
    fs["adaptiveThreshWinSizeMin"] >> params->adaptiveThreshWinSizeMin;
    fs["adaptiveThreshWinSizeMax"] >> params->adaptiveThreshWinSizeMax;
    fs["adaptiveThreshWinSizeStep"] >> params->adaptiveThreshWinSizeStep;
    fs["adaptiveThreshConstant"] >> params->adaptiveThreshConstant;
    fs["minMarkerPerimeterRate"] >> params->minMarkerPerimeterRate;
    fs["maxMarkerPerimeterRate"] >> params->maxMarkerPerimeterRate;
    fs["polygonalApproxAccuracyRate"] >> params->polygonalApproxAccuracyRate;
    fs["minCornerDistanceRate"] >> params->minCornerDistanceRate;
    fs["minDistanceToBorder"] >> params->minDistanceToBorder;
    fs["minMarkerDistanceRate"] >> params->minMarkerDistanceRate;
    //fs["cornerRefinementMethod"] >> params->cornerRefinementMethod;
    fs["cornerRefinementWinSize"] >> params->cornerRefinementWinSize;
    fs["cornerRefinementMaxIterations"] >> params->cornerRefinementMaxIterations;
    fs["cornerRefinementMinAccuracy"] >> params->cornerRefinementMinAccuracy;
    fs["markerBorderBits"] >> params->markerBorderBits;
    fs["perspectiveRemovePixelPerCell"] >> params->perspectiveRemovePixelPerCell;
    fs["perspectiveRemoveIgnoredMarginPerCell"] >> params->perspectiveRemoveIgnoredMarginPerCell;
    fs["maxErroneousBitsInBorderRate"] >> params->maxErroneousBitsInBorderRate;
    fs["minOtsuStdDev"] >> params->minOtsuStdDev;
    fs["errorCorrectionRate"] >> params->errorCorrectionRate;
    return true;
}


/*
 * Helper function that returns the sign of a float.
 */
static inline int sign(float f)
{
    if (f > 0) return 1;
    if (f < 0) return -1;
    return 0;
}


/*
 * Squared length for a point / vector.
 */
static inline float distsq(const Point2f &p)
{
    return p.x * p.x + p.y * p.y;
}


// Point in triangle test, adapted from:
// https://stackoverflow.com/questions/2049582/how-to-determine-if-a-point-is-in-a-2d-triangle
float side (const Point2f &p1, const Point2f &p2, const Point2f &p3)
{
    return (p1.x - p3.x) * (p2.y - p3.y) - (p2.x - p3.x) * (p1.y - p3.y);
}

bool pointInTriangle (const Point2f &pt, const Point2f &v1, const Point2f &v2, const Point2f &v3)
{
    bool b1, b2, b3;

    b1 = side(pt, v1, v2) < 0.0f;
    b2 = side(pt, v2, v3) < 0.0f;
    b3 = side(pt, v3, v1) < 0.0f;

    return ((b1 == b2) && (b2 == b3));
}


/*
 * Helper function to check if two floats are (nearly) equal.
 */
bool flt_eq(float a, float b)
{
    const float max_diff = 0.00001;
    float d = a - b;
    if (d < 0) d = -d;
    return d < max_diff;
}


/*
 * Helper function to transform a point with a matrix.
 */
static Point2f transformPoint(const Point2f &p, const Mat &t)
{
    Mat m(3, 1, CV_64F);
    m.at<double>(0, 0) = p.x;
    m.at<double>(1, 0) = p.y;
    m.at<double>(2, 0) = 1;

    Mat r = t * m;

    double w = r.at<double>(0, 2);
    return Point2f(r.at<double>(0, 0) / w, r.at<double>(0, 1) / w);
}


/*
 * This generates a vector of quads for a given set of image and world points.
 * Matching points in p_img and p_world must have the same index.
 * grid_size specifies the length of square edge in worldspace.
 */
static vector<MappedQuad> calculateQuads(const vector<Point2f> &p_img, const vector<Point2f> &p_world, float grid_size)
{
    vector<MappedQuad> ret;

    for (size_t i = 0; i < p_world.size(); i++)
    {
        // For every point: Try to find up to find the "forward quad" (The quad with all other points > current)
        size_t points[4];
        points[0] = i;
        points[1] = SIZE_MAX;
        points[2] = SIZE_MAX;
        points[3] = SIZE_MAX;

        int cnt = 0;

        // The position of the current point
        float xa = p_world[i].x;
        float ya = p_world[i].y;

        for (size_t j = 0; j < p_world.size(); j++)
        {
            // The position of the other point
            float xb = p_world[j].x;
            float yb = p_world[j].y;

            if (flt_eq(xa + grid_size, xb))
            {
                // Point is in next column -> check row
                if (flt_eq(ya + grid_size, yb))
                {
                    // -> row below
                    points[3] = j;
                    cnt++;
                }
                else if (flt_eq(ya, yb))
                {
                    // -> current row
                    points[1] = j;
                    cnt++;
                }
            }
            else if (flt_eq(ya + grid_size, yb) && flt_eq(xa, xb))
            {
                // current column and row below
                points[2] = j;
                cnt++;
            }

            if (cnt == 3)
            {
                // No need to continue when I have all 3 points
                break;
            }
        }

        if (cnt != 3)
        {
            //cout << "No quad for " << xa << ", " << ya << endl;
            // No Quad found
            continue;
        }

        // OK, got a quad -> get the points and calculate transformation
        MappedQuad m;
        Point2f wld[4];

        for (size_t j = 0; j < 4; j++)
        {
            if (points[j] == SIZE_MAX)
            {
                // this should never happen
                cout << "DUPLICATE POINTS!" << endl;
                cerr << "DUPLICATE POINTS!" << endl;
                exit(1);
                return ret;
            }

            // copy image points to quad
            m.corners_image[j] = p_img[points[j]];

            wld[j] = p_world[points[j]];

            //cout << m.corners_image[j].x << ", " << m.corners_image[j].y << " -> " << wld[j].x << ", " << wld[j].y << endl;
        }

        // Edge points
        m.ws_min = wld[0];
        m.ws_max = wld[3];

        // Initialize adjacent_quads with invalid value
        for (size_t j = 0; j < 8; j++)
        {
            m.adjacent_quads[j] = SIZE_MAX;
        }

        // calculate transform
        m.transform = cv::getPerspectiveTransform(m.corners_image, wld);

        /*
        for (size_t j = 0; j < 4; j++)
        {
            Point2f r = transformPoint(m.corners_image[j], m.transform);
            cout << "RES: " << r.x << ", " << r.y << endl;
        }
        */

        ret.push_back(m);
    }
    return ret;
}


/*
 * Finds the adjacent quads for all quads in a given set of quads.
 */
void findAdjacentQuads(vector<MappedQuad> &quads)
{

    /*
     * Edges of a quad:
     *
     * 0==1
     * |  |
     * 2==3
     *
     * Quads:
     *
     *  *---*---*---*
     *  : 0 : 1 : 2 :
     *  *---*===*---*
     *  : 3 |(8)| 4 :
     *  *---*===*---*
     *  : 5 : 6 : 7 :
     *  *---*---*---*
     *
     * Assuming that the input is sane, I only need one corner compare to detect a neighbor.
     */

    for (size_t i = 0; i < quads.size(); i++)
    {
        auto &q = quads[i];

        // Find neighbors for q
        for (size_t j = 0; j < quads.size(); j++)
        {
            const auto &w = quads[j];
            if      (q.corners_image[0] == w.corners_image[3]) q.adjacent_quads[0] = j;
            else if (q.corners_image[0] == w.corners_image[1]) q.adjacent_quads[3] = j;
            else if (q.corners_image[0] == w.corners_image[2]) q.adjacent_quads[1] = j;
            else if (q.corners_image[1] == w.corners_image[2]) q.adjacent_quads[2] = j;
            else if (q.corners_image[1] == w.corners_image[0]) q.adjacent_quads[4] = j;
            else if (q.corners_image[2] == w.corners_image[1]) q.adjacent_quads[5] = j;
            else if (q.corners_image[2] == w.corners_image[0]) q.adjacent_quads[6] = j;
            else if (q.corners_image[3] == w.corners_image[0]) q.adjacent_quads[7] = j;
        }
    }
}


/*
 * Finds the quad containing a given point.
 * The <hint> is the index where the search starts. If the last found quad is given as hint and <p> is next
 * to the last point, the search will immediately find the quad containing <p>.
 */
static size_t findQuad(const vector<MappedQuad> &quads, const Point2f &p, size_t hint)
{
    if (quads.empty())
    {
        return SIZE_MAX;
    }

    size_t i = hint;

    if (i > quads.size())
    {
        i = 0;
    }

    // Loop over all quads and make a point-in-triangle check for all image triangles.
    do
    {
        const auto &q = quads[i];
        if (pointInTriangle(p, q.corners_image[0], q.corners_image[1], q.corners_image[2]) ||
            pointInTriangle(p, q.corners_image[3], q.corners_image[1], q.corners_image[2]))
        {
            return i;
        }
        i++;
        if (i == quads.size())
        {
            i = 0;
        }

    } while (i != hint);
    return SIZE_MAX;
}


/*
 * Calculates the overlap factor for a quad form (point - size) to (point + size) and the quad q.
 * The overlap factor is the portion of area in <q> covered by the quad around <point>.
 * It is always between 0 (no overlap) and 1 (full overlap).
 */
static float calculateWsOverlap(const MappedQuad &q, const Point2f &point, float size)
{
    float sqs = q.ws_max.x - q.ws_min.x;
    size *= sqs;
    //cout << "overlap for (" << point.x << "," << point.y << ") min=(" << q.ws_min.x << ", " << q.ws_min.y << ") max=(" << q.ws_max.x << ", " << q.ws_max.y << ")" << endl;
    float xl = point.x - size;
    float xh = point.x + size;

    if (xl > q.ws_max.x || xh < q.ws_min.x)
        return 0;

    // Clamp to bounds
    if (xl < q.ws_min.x) xl = q.ws_min.x;
    if (xh > q.ws_max.x) xh = q.ws_max.x;
    //cout << "xl=" << xl << " xh=" << xh << endl;

    float yl = point.y - size;
    float yh = point.y + size;

    if (yl > q.ws_max.y || yh < q.ws_min.y)
        return 0;

    // Clamp to bounds
    if (yl < q.ws_min.y) yl = q.ws_min.y;
    if (yh > q.ws_max.y) yh = q.ws_max.y;
    //cout << "yl=" << yl << " yh=" << yh << endl;

    float w = xh - xl;
    float h = yh - yl;

    //cout << "w=" << w << " h=" << h << endl;

    if (w <= 0 || h <= 0)
        return 0;

    return (w / sqs) * (h / sqs);
}

/*
 * Find the <n> quads that are closest to a given image space point.
 * max_dist limit the maximal distance for a quad.
 * (The average of the squared distances form all points must be smaller than max_dist)
 *
 * The return value is a map form distance to the point index.
 */
static map<float, size_t> findNearestQuads(const vector<MappedQuad> &quads, const Point2f &center, int n, int max_dist)
{
    // I'm using a map here as a sorted container for all the distances
    map<float, size_t> min_quads;

    for (size_t i = 0; i < quads.size(); i++) {
        const auto &q = quads[i];
        float dist =
            distsq(q.corners_image[0] - center) +
            distsq(q.corners_image[1] - center) +
            distsq(q.corners_image[2] - center) +
            distsq(q.corners_image[3] - center);

        if (dist > max_dist * 4)
        {
            // reject, too far away
            continue;
        }

        // Check if empty or last element is larger than current.
        if (min_quads.empty() || (min_quads.crbegin()->first > dist)) {
            min_quads[dist] = i;

            // if i have more than n elements delete the last
            if (min_quads.size() > n) {
                min_quads.erase(prev(min_quads.end()));
            }
        }
    }

    return min_quads;
}


/*
 * Calculates a worldspace point for a given images point by extrapolating from the nearest quads.
 *
 * z           z-coordinate for this point.
 * max_dist    maximum distance for the quads used for extrapolation
 * min_quads   minimum number of quads in range to do the extrapolation
 * max_err     maximum difference between the the extrapolation results for all quads
 *
 * If the extrapolation fails, the w coordinate of the returned point is -1
 */
static Vec4f extraploatePoint(const vector<MappedQuad> &quads, const Point2f &center, float z, int max_dist, int min_quads, float max_err)
{

    /*
     * The minimum quality value of an extrapolated pixel.
     */
    const float MIN_QUALITY = 1000.0f;

    /*
     * Number of quads used for extrapolation
     */
    const int exp_max = std::max(5, min_quads * 2);

    // Find the quads
    const auto source = findNearestQuads(quads, center, exp_max, max_dist);

    // check if enough
    if (source.empty() || source.size() < min_quads)
    {
        //cout << "Extrapolation failed! no mapped quad!" << endl;
        return Vec4f(0, 0, 0, -1);
    }

    //cout << "Doing extrapolation with " << source.size() << " Points" << endl;

    // For all source quads transform the pixel to the world space and then calculate the average by using the distances as weighting factor.
    Point2f sum_vect = Vec2f(0, 0);
    float sum = 0;
    float quality = 0;

    for (const auto &s : source) {
        const auto &a = quads[s.second];
        Point2f w = transformPoint(center, a.transform);

        if (sum != 0)
        {
            // The pixel is NOT mapped if the extrapolation results for the quads are too far from each other
            Point2f avg = sum_vect * (1 / sum);
            float ad = distsq(avg - w);
            if (ad > (max_err * max_err))
            {
                //cout << "Discard point with distance of " << ad << " to the avgerage as garbage" << endl;
                return Vec4f(0, 0, 0, -1);
            }
        }

	    float f = 1 / s.first;
        sum_vect += w * f;
        sum += f;
        quality += s.first;
    }

    // Use the disatnce as quality value for the pixel.
    // By doing so I ensure, that an extrapolated point has not too much influence on a good ray.
    quality /= source.size();
    if (quality < MIN_QUALITY)
    {
        quality = MIN_QUALITY;
    }

    return Vec4f(sum_vect.x / sum, sum_vect.y / sum, z, quality);
}

/*
 * Interpolates a point for a given set of mapped quads
 *
 * z           z-coordinate for this point
 * hint        the quad containing the last point. this is set to the current quad when a quad is found
 *
 * max_exp_dist, min_exp_quads, max_exp_err
 * Parameters for the extrapolation
 *
 * The z coordinate is set to the given value.
 * the w coordinate of the returned point is the "quality" of the returned point.
 * If the interpolation fails, the w coordinbate is negative
 */
static Vec4f interpolatePoint(const vector<MappedQuad> &quads, const Point2f &center, float z, size_t &hint, int max_exp_dist, int min_exp_quads, float max_exp_err)
{
    /*
     * Find the quad containing the point.
     * Then use the transform to project the point into the worldspace.
     * Using this worldspace positions, the weighting factors for the adjacent quads are
     * calculated. The final worldspace pos is the weighted average of the nine quads.
     */



    /*
     * This is the size in each direction from the center of the interpolation rectangle.
     * IP_SIZE * 2 is the edge length of the total area used for interpolation.
     * It should must be smaller or equal to 1 (1 full quad size in each direction).
     * If the value is < 0.5 only the current quad is used for the center pixels.
     */
    const float IP_SIZE = 1.4f / 2;

    // Find the current quad
    size_t quad_idx = findQuad(quads, center, hint);

    if (quad_idx == SIZE_MAX)
    {
        // No quad found -> extrapolate
        if (max_exp_dist <= 0)
        {
            // No extrapolation
            return Vec4f(0, 0, 0, -1);
        }
        return extraploatePoint(quads, center, z, max_exp_dist, min_exp_quads, max_exp_err);
    }

    // update the hint
    hint = quad_idx;

    //cout << center.x << ", " << center.y << " is on " << quad_idx << endl;

    const auto &q = quads[quad_idx];

    // ws is the point in worldspace
    Point2f ws = transformPoint(center, q.transform);

    // Sum for averaging
    Point2f sum_vect = Vec2f(0, 0);
    float sum_div = 0;

    // Minimum (squared) distance between the point and a edge
    float min_dist = FLT_MAX;

    for (int i = 0; i < 9; i++)
    {
        if (i != 8 && q.adjacent_quads[i] == SIZE_MAX)
        {
            // No neighbour at this position
            continue;
        }
        const auto &a = (i == 8) ? q : quads[q.adjacent_quads[i]];

        float factor = calculateWsOverlap(a, ws, IP_SIZE);

        if (factor <= 0)
        {
            // Quad is not used (no overlap)
            continue;
        }

        // Calculate the minimum distance to a vertex.
        // This is used as a quality indicator
        for (int i = 0; i < 4; i++)
        {
            float d = distsq(a.corners_image[i] - center);
            if (d < min_dist)
            {
                min_dist = d;
            }
        }

        Point2f w = transformPoint(center, a.transform);
        //cout << center.x << ", " << center.y << " -> " << i << ": " << w.x << ", " << w.y << " / " << factor << endl;
        sum_vect += w * factor;
        sum_div += factor;
    }

    return Vec4f(sum_vect.x / sum_div, sum_vect.y / sum_div, z, std::sqrt(min_dist));
}


/*
 * Calculate a ray starting at the 0 plane form two points.
 */
static EstimatedRay calculateNormalizedRay(const Vec4f &a, const Vec4f &b)
{
    // Get rid of the w coordinate and calculate distance
    Vec3f start = Vec3f(a[0], a[1], a[2]);
    Vec3f diff(b[0] - a[0], b[1] - a[1], b[2] - a[2]);

    // Normalize distance
    diff = diff / norm(diff);

    // Project start point onto 0 plane
    float dist = start[2] / diff[2];
    start -= diff * dist;


    EstimatedRay er;
    er.start = start;
    er.dir = diff;

    // Better quality if distance is greater
    er.rating = (a[3] + b[3]) / (b[2] - a[2]);

    if (er.start[0] != er.start[0] ||
        er.start[0] != er.start[0] ||
        er.start[0] != er.start[0])
    {
        cout << "a=" << a << " b=" << b << " diff=" << diff << " dist=" << dist << endl;
        exit(0);
    }

    return er;
}

/*
 * Calculate the point where two lines are nearest
 * see: https://math.stackexchange.com/questions/1036959/midpoint-of-the-shortest-distance-between-2-rays-in-3d?utm_medium=organic&utm_source=google_rich_qa&utm_campaign=google_rich_qa
 */
static Point3f estimateConvergencePoint(const EstimatedRay &ra, const EstimatedRay &rb)
{
    Point3f r1o(ra.start[0], ra.start[1], ra.start[2]);
    Point3f r2o(rb.start[0], rb.start[1], rb.start[2]);
    Point3f r1d(ra.dir[0], ra.dir[1], ra.dir[2]);
    Point3f r2d(rb.dir[0], rb.dir[1], rb.dir[2]);

    float t =
        (((r2o - r1o).dot(r1d)) * r2d.dot(r2d) + (r1o - r2o).dot(r2d) * r1d.dot(r2d)) /
        (r1d.dot(r1d) * r2d.dot(r2d) - (r1d.dot(r2d) * r1d.dot(r2d)));
    float s =
        (((r1o - r2o).dot(r2d)) * r1d.dot(r1d) + (r2o - r1o).dot(r1d) * r1d.dot(r2d)) /
        (r1d.dot(r1d) * r2d.dot(r2d) - (r1d.dot(r2d) * r1d.dot(r2d)));

    Point3f m = ((r1o + r1d * t) + (r2o + r2d * s)) / 2;
    return m;
}


/*
 * Calculates the average point of a set of points.
 * The worst (farthest form the average) <reject> points are deleted from the list before calculating the final average.
 */
static Point3f calculateAveragePoint(vector<Point3f> &points, size_t reject, float &err)
{
    Point3f avg;
    if (reject != 0)
    {
        if (reject > points.size())
        {
            err = 0;
            return Point3f();
        }

        // Calculate the initial average
        for (const auto &p : points)
        {
            avg += p;
        }
        avg *= 1.0f / points.size();


        vector<pair<float, size_t>> v;
        v.reserve(points.size());

        // Sort points by distance to average
        for (size_t i = 0; i < points.size(); i++)
        {
            float d = norm(avg - points[i]);
            v.push_back(make_pair(d, i));
        }
        sort(v.begin(), v.end());

        // Copy the (size - reject) best points to res
        vector<Point3f> res;
        res.resize(v.size() - reject);
        for (size_t i = 0; i < res.size(); i++)
        {
            res[i] = points[v[i].second];
        }
        cout << "num of points: " << points.size() << " worst: " << v[v.size() - 1].first  << "idx: " << v[v.size() - 1].second << " last used: " << v[res.size() - 1].first << " idx: " << v[res.size() - 1].second << endl;

        // overwrite old point list
        points = res;
    }

    // Now 'points' only contains the better points -> recalculate average and error
    avg = Point3f();
    err = 0;

    for (const auto &p : points)
    {
        avg += p;
    }

    avg *= 1.0f / points.size();

    for (const auto &p : points)
    {
        float d = norm(avg - p);
        err += d;
    }
    err /= points.size();

    return avg;
}

/**
 */
int main(int argc, char *argv[]) {
    /*
    EstimatedRay ra;
    ra.start = Vec3f(0.25485, 0.15796, 0);
    ra.dir = Vec3f(-0.0446953, -0.118038, 0.992003);

    EstimatedRay rb;
    rb.start = Vec3f(.529294, 0.590926, 0);
    rb.dir = Vec3f(-0.326208, -0.554728, 0.765418);

    cout << estimateConvergencePoint(ra, rb) << endl;

    rb.start += rb.dir * 100;
    cout << estimateConvergencePoint(ra, rb) << endl;

    return 0;
    */

    CommandLineParser parser(argc, argv, keys);
    parser.about(about);

    if(argc < 6) {
        parser.printMessage();
        return 0;
    }

    int squaresX = parser.get<int>("w");
    int squaresY = parser.get<int>("h");
    float squareLength = parser.get<float>("sl");
    float markerLength = parser.get<float>("ml");
    int dictionaryId = parser.get<int>("d");
    bool refindStrategy = parser.has("rs");
    String input = parser.get<String>("v");
    String output = parser.get<String>("o");

    int extrapolation_dist = 200;
    if (parser.has("ed"))
    {
        extrapolation_dist = parser.get<int>("ed");
    }

    int extrapolation_min_quads = parser.get<int>("em");
    float extrapolation_max_err = parser.get<float>("ee");

    // This is compared to squared distance values so square it.
    extrapolation_dist *= extrapolation_dist;

    float dist = 1.0f;
    if (parser.has("z"))
    {
        dist = parser.get<float>("z");
    }

    cout << "ed=" << extrapolation_dist << " em=" << extrapolation_min_quads << " ee=" << extrapolation_max_err << endl;

    Ptr<aruco::DetectorParameters> detectorParams = aruco::DetectorParameters::create();
	detectorParams->adaptiveThreshWinSizeMin = 3; //3
    detectorParams->adaptiveThreshWinSizeMax = 30; //23
    detectorParams->adaptiveThreshWinSizeStep = 5; //10

    //detectorParams->adaptiveThreshConstant = 7; //7
/*    params->minMarkerPerimeterRate;
    params->maxMarkerPerimeterRate;
    params->polygonalApproxAccuracyRate;
    params->minCornerDistanceRate;
    params->minDistanceToBorder;
    params->minMarkerDistanceRate;
    params->cornerRefinementWinSize;
    params->cornerRefinementMaxIterations;
    params->cornerRefinementMinAccuracy;
    params->markerBorderBits;
    params->perspectiveRemovePixelPerCell;
    params->perspectiveRemoveIgnoredMarginPerCell;
    params->maxErroneousBitsInBorderRate;
    params->minOtsuStdDev;
    params->errorCorrectionRate;*/
	
	//cout << "after readDetectorParameters "<<detectorParams->adaptiveThreshWinSizeMin<<" "<<detectorParams->adaptiveThreshWinSizeMax<<" "<<detectorParams->adaptiveThreshWinSizeStep<<" "<<detectorParams->adaptiveThreshConstant<<"\n";
	
	
	
	

    if(parser.has("dp")) {
		cout << "in readDetectorParameters "<<parser.get<string>("dp")<<" \n";
        bool readOk = readDetectorParameters(parser.get<string>("dp"), detectorParams);
        if(!readOk) {
            cerr << "Invalid detector parameters file" << endl;
            return 0;
        }
    }

    if(!parser.check()) {
        parser.printErrors();
        return 0;
    }

    Ptr<aruco::Dictionary> dictionary =
        aruco::getPredefinedDictionary(aruco::PREDEFINED_DICTIONARY_NAME(dictionaryId));

    VideoCapture inputVideo;
    inputVideo.open(input);

    float axisLength = 0.5f * ((float)min(squaresX, squaresY) * (squareLength));

    // create charuco board object
    Ptr<aruco::CharucoBoard> charucoboard =
        aruco::CharucoBoard::create(squaresX, squaresY, squareLength, markerLength, dictionary);
    Ptr<aruco::Board> board = charucoboard.staticCast<aruco::Board>();

    double totalTime = 0;
    int totalIterations = 0;

    //namedWindow("out",WINDOW_NORMAL);

    Mat camMatrix, distCoeffs;

    typedef vector<vector<Vec4f>> PixelMap;

    // This will contain all worldspace points for all image points (in all images)
    vector<PixelMap> pixel_to_world;

    // The current z value
    float z = 0;

    Mat firstImage;

    while(inputVideo.grab()) {
        Mat image, imageCopy;
        inputVideo.retrieve(image);
    	

        if (firstImage.size().width == 0)
        {
            firstImage = image;
        }

        cout << "process image at z=" << z << endl;

        double tick = (double)getTickCount();

        vector< int > markerIds, charucoIds;
        vector< vector< Point2f > > markerCorners, rejectedMarkers;
        vector< Point2f > charucoCorners;
        Vec3d rvec, tvec;

        // detect markers
        aruco::detectMarkers(image, dictionary, markerCorners, markerIds, detectorParams,
                             rejectedMarkers);

        // refind strategy to detect more markers
        //if(refindStrategy)
		aruco::refineDetectedMarkers(image, board, markerCorners, markerIds, rejectedMarkers, camMatrix, distCoeffs);
//Ptr<aruco::Board> board2 = aruco::GridBoard::create(squaresX, squaresY, squareLength, markerLength, dictionary);
//aruco::refineDetectedMarkers(image, board2, markerCorners, markerIds, rejectedMarkers);

        // interpolate charuco corners
        int interpolatedCorners = 0;
        if(markerIds.size() > 0)
            interpolatedCorners =
                aruco::interpolateCornersCharuco(markerCorners, markerIds, image, charucoboard,
                                                 charucoCorners, charucoIds, camMatrix, distCoeffs);



        double currentTime = ((double)getTickCount() - tick) / getTickFrequency();
        totalTime += currentTime;
        totalIterations++;
        if(totalIterations % 30 == 0) {
            cout << "Detection Time = " << currentTime * 1000 << " ms "
                 << "(Mean = " << 1000 * totalTime / double(totalIterations) << " ms)" << endl;
        }



		cout << "markers size = " << markerIds.size() << ", "<< charucoIds.size()<<", "<<interpolatedCorners << endl;
		
// id erkannt? korrekte interpretation?; print IDs
// aruco::drawDetectedCornersCharuco(imageCopy, currentCharucoCorners, currentCharucoIds);
// interpolateCornersCharuco   -->~board?

// points the same?

	
        // draw results
        image.copyTo(imageCopy);
        if(markerIds.size() > 0) {
            aruco::drawDetectedMarkers(imageCopy, markerCorners);
			char buff2[1000];
			snprintf(buff2, sizeof(buff2), (output.operator std::string() + "_markers_%1.5f.png").c_str(), z);	
			imwrite(buff2, imageCopy);        
        }

        if(interpolatedCorners > 0) {
            Scalar color;
            color = Scalar(255, 0, 0);
            aruco::drawDetectedCornersCharuco(imageCopy, charucoCorners, charucoIds, color);
			char buff3[1000];
			snprintf(buff3, sizeof(buff3), (output.operator std::string() + "_markersC_%1.5f.png").c_str(), z);	
			imwrite(buff3, imageCopy);
        
        }

		//imshow("out", imageCopy);
        //while(waitKey(0) != 32);

        if (charucoIds.size() != charucoCorners.size())
        {
            cout << "Length of id list must be equal to length of corner list" << endl;
            return 1;
        }

        vector<Point2f> dataPoints;
        dataPoints.resize(charucoIds.size());

        for (size_t corner = 0; corner < charucoIds.size(); corner++)
        {
            const auto &w = charucoboard->chessboardCorners[charucoIds[corner]];
            //cout << "W=[" << w.x << ", " << w.y << ", " << z << "] -> P=";
            //cout << charucoCorners[corner] << endl;
            dataPoints[corner] = Point2f(w.x, w.y);
        }

        // Add new entry (with the size of the image) to pixel mapping
        pixel_to_world.push_back(PixelMap(image.size().height, vector<Vec4f>(image.size().width, Vec4f())));
        auto &interpolated_pix_pos = pixel_to_world[pixel_to_world.size() - 1];


        // Find the quads
        cout << "Make quads" << endl;
        vector<MappedQuad> quads = calculateQuads(charucoCorners, dataPoints, squareLength);
        cout << "Found " << quads.size() << " quads" << endl;
        findAdjacentQuads(quads);

        // Run interpolation (and extrapolation) for all pixels
        // This is done in multiple threads
        cout << "Interpolate points" << endl;
        size_t height = interpolated_pix_pos.size();
        tbb::parallel_for(size_t(0), height, size_t(1), [=,&interpolated_pix_pos](size_t y)
        {
            size_t hint = 0;
            auto &row = interpolated_pix_pos[y];
            for (size_t x = 0; x < row.size(); x++)
            {
                row[x] = interpolatePoint(quads, Point2f(x, y), z, hint, extrapolation_dist, extrapolation_min_quads, extrapolation_max_err);

                /*
                if (x % 50 == 0 && y % 50 == 0 && row[x][3] > 0)
                {
                    cout << "(" << x << ", " << y << ") -> (" << row[x][0] << ", " << row[x][1] << ", " << row[x][2] << ")" << endl;
                }
                */
            }
        });

        // Debug images for the single steps
        // These image simply encode the worldspace positions in the colors.
        // For an unmapped pixel the input is copied.
        #if 1
        Mat distImg(interpolated_pix_pos.size(), interpolated_pix_pos[0].size(), CV_8UC3);
        for (size_t y = 0; y < interpolated_pix_pos.size(); y++)
        {
            auto &row = interpolated_pix_pos[y];
            for (size_t x = 0; x < row.size(); x++)
            {
                Vec3b &d = distImg.at<Vec3b>(y, x);
                const Vec3b &s = image.at<Vec3b>(y, x);

                if (row[x][3] < 0)
                {
                    d = s;
                }
                else
                {
                    float dx = 0;
                    if (x < row.size() - 1)
                    {
                        dx = row[x][0];
                    }
                    float dy = 0;
                    if (y < interpolated_pix_pos.size() - 1)
                    {
                        dy = row[x][1];
                    }

                    dx *= 500;
                    dy *= 500;

                    if (dx < 0)
                    {
                        dx = 0;
                    }

                    if (dy < 0)
                    {
                        dy = 0;
                    }

                    if (dx > 255)
                    {
                        dx = 255;
                    }

                    if (dy > 255)
                    {
                        dy = 255;
                    }

                    d[0] = dx;
                    d[1] = dy;
                    d[2] = s[0];
                }
            }
        }

		char buff[1000];
		snprintf(buff, sizeof(buff), (output.operator std::string() + "_dist_%1.5f.png").c_str(), z);	
        imwrite(buff, distImg);
        //cout << "show dist img, x=" << distImg.size().width << " inp.x=" << image.size().width << endl;
        //imshow("out", distImg);
        //while(waitKey(0) != 32);
        #endif

        cout << "Next frame" << endl;
        z += dist;
    }

    /*
     * Now we have worldspace positions for all mapped pixels.
     * Its time to calculate the rays form all these pixels.
     */

    cout << "Done processing images -> calculate rays" << endl;


    // This will contain the final rays
    vector<vector<EstimatedRay>> final_rays(pixel_to_world[0].size(), vector<EstimatedRay>(pixel_to_world[0][0].size(), EstimatedRay()));

    // Error values for debugging and to check the results
    float max_err_start = 0;
    float max_err_dir = 0;
    float avg_err_start = 0;
    float avg_err_dir = 0;

    // Weighted errors
    float wavg_err_start = 0;
    float wavg_err_dir = 0;

    int total = 0;

    // Min / max direction vectors
    float min_dx = FLT_MAX, max_dx = FLT_MIN;
    float min_dy = FLT_MAX, max_dy = FLT_MIN;

    // The error image shows the differences between the sub-rays for all pixels.
    Mat err_img(final_rays.size(), final_rays[0].size(), CV_32FC3);

    for (size_t y = 0; y < final_rays.size(); y++)
    {
        for (size_t x = 0; x < final_rays[y].size(); x++)
        {
            // All rays between all images for this pixel.
            
            vector<EstimatedRay> rays;

            for (size_t a = 0; a < pixel_to_world.size(); a++)
            {
				bool aIsNan = pixel_to_world[a][y][x][0] != pixel_to_world[a][y][x][0] || 
						pixel_to_world[a][y][x][1] != pixel_to_world[a][y][x][1] || 
						pixel_to_world[a][y][x][2] != pixel_to_world[a][y][x][2];
				if(aIsNan){
					//cout << "encountered nan: a=" << pixel_to_world[a][y][x] << endl;
				}
				else{
					for (size_t b = a + 1; b < pixel_to_world.size(); b++)
					{
						bool bIsNan = pixel_to_world[b][y][x][0] != pixel_to_world[b][y][x][0] || 
							pixel_to_world[b][y][x][1] != pixel_to_world[b][y][x][1] || 
							pixel_to_world[b][y][x][2] != pixel_to_world[b][y][x][2];
						// Calculate rays from a to b plane
						if (pixel_to_world[a][y][x][3] < 0 || pixel_to_world[b][y][x][3] < 0 || bIsNan)
						{
							// a or b are not mapped
							if(bIsNan){
								//cout << "encountered nan: b=" << pixel_to_world[b][y][x] << endl;
							}
						}else{
							EstimatedRay er = calculateNormalizedRay(pixel_to_world[a][y][x], pixel_to_world[b][y][x]);
							rays.push_back(er);
						}
					}
				}
            }

            // Now 'rays' contains all rays in a normalized form.
            // The final rays are calculated by averaging over the rays

            // Set error image to black at this point
            Vec3f &d = err_img.at<Vec3f>(y, x);
            d = Vec3f(0, 0, 0);

            auto &final_ray = final_rays[y][x];
            if (rays.size() == 0)
            {
                // Not a single ray for this pixel -> mark as unmapped
                final_ray.rating = -1;
                continue;
            }

            // Sort the rays by the rating
            sort(rays.begin(), rays.end());

            if (rays.size() > 50)
            {
                // only use the best 50 rays
                rays.resize(50);
            }
            // Maybe also discard the worst rays if they are extraploated and the best is not?

			// add check for excluding really bad rays! overhead but should work
			{
				float total_rating = 0;
				for (const auto &r : rays)
				{
					total_rating += 1 / r.rating;
				}
				
				EstimatedRay final_ray_dummy;
				// Weighted sum of the sub-rays.
				for (const auto &r : rays)
				{
					final_ray_dummy.start += (r.start * (1 / r.rating)) / total_rating;
					final_ray_dummy.dir += (r.dir * (1 / r.rating));
				}

				// normalize the direction vector
				final_ray_dummy.dir /= norm(final_ray_dummy.dir);
				
				int rayIndex = 0;
				int raysSize = rays.size();
				for (int i = 0; i < raysSize; i++)
				{
					float es = norm(final_ray_dummy.start - rays[rayIndex].start);
					float ed = norm(final_ray_dummy.dir - rays[rayIndex].dir);

					if(es>5 || ed>0.5){
						rays.erase(rays.begin()+rayIndex);
					}else{
						rayIndex++;
					}
				}
				if (rays.size() == 0){
					cout << "zero rays left - should not happen! what to do?" <<endl;
					final_ray.rating = -1;
					final_ray.start = Vec3f(0,0,0);
					final_ray.dir   = Vec3f(0,0,0);
					d = Vec3f(0, 0, 0);
					total--;
					continue;
				}
			}


            // Mapped -> set blue to 1 (max)
            d[0] = 1;

            total++;

            // Calculate the total rating as sum of all 1 / rating.
            float total_rating = 0;
            for (const auto &r : rays)
            {
                total_rating += 1 / r.rating;
            }

            // Weighted sum of the sub-rays.
            for (const auto &r : rays)
            {
                final_ray.start += (r.start * (1 / r.rating)) / total_rating;
                final_ray.dir += (r.dir * (1 / r.rating));
            }

            // normalize the direction vector
            final_ray.dir /= norm(final_ray.dir);

            // Calculate errors
            float es_sum = 0, ed_sum = 0;
            float es_wsum = 0, ed_wsum = 0;
            for (const auto &r : rays)
            {
                float es = norm(final_ray.start - r.start);
                float ed = norm(final_ray.dir - r.dir);

                if (es > max_err_start) max_err_start = es;
                if (ed > max_err_dir) max_err_dir = ed;

                es_sum += es;
                ed_sum += ed;

                es_wsum += es / r.rating;
                ed_wsum += ed / r.rating;
            }

            //cout << "num of rays: " << rays.size() << "rating: " << total_rating << "frs: " << final_ray.start << "r0s" << rays[0].start << endl;

            es_sum /= rays.size();
            ed_sum /= rays.size();

            es_wsum /= total_rating;
            ed_wsum /= total_rating;

            avg_err_start += es_sum;
            avg_err_dir += ed_sum;

            wavg_err_start += es_wsum;
            wavg_err_dir += ed_wsum;

            d[1] = es_sum;
            d[2] = ed_sum;

            if (final_ray.dir[0] > max_dx) max_dx = final_ray.dir[0];
            if (final_ray.dir[0] < min_dx) min_dx = final_ray.dir[0];
            if (final_ray.dir[1] > max_dy) max_dy = final_ray.dir[1];
            if (final_ray.dir[1] < min_dy) min_dy = final_ray.dir[1];

            final_ray.rating = total_rating;
        }
    }

    {
        // normalize error image
        Mat normalized_err_img(err_img.size().height, err_img.size().width, CV_8UC3);
        for (size_t y = 0; y < err_img.size().height; y++)
        {
            for (size_t x = 0; x < err_img.size().width; x++)
            {
                const Vec3f &d = err_img.at<Vec3f>(y, x);
                Vec3b &n = normalized_err_img.at<Vec3b>(y, x);
                n[0] = 0;
                n[1] = (d[1] / max_err_start) * 255;
                n[2] = (d[2] / max_err_dir) * 255;
            }
        }

        imwrite(output + "_err.png", normalized_err_img);
        //imshow("out", normalized_err_img);
        //while(waitKey(0) != 32);
    }

    if (total != 0)
    {
        avg_err_start /= total;
        avg_err_dir /= total;
        wavg_err_start /= total;
        wavg_err_dir /= total;
    }

    ostringstream errorstr;

    errorstr << "Errors: max_start=" << max_err_start << " max_dir=" << max_err_dir << " avg_start=" << avg_err_start << " avg_dir=" << avg_err_dir << endl;
    errorstr << "        wavg_start=" << wavg_err_start << " wavg_dir=" << wavg_err_dir << endl;

    Mat dir_img(final_rays.size(), final_rays[0].size(), CV_8UC3);

    // OK, writing output file
    ofstream result;
    result.open (output + "_rays.txt", ios::out | ios::trunc);

    result << final_rays.size() << ' ' << final_rays[0].size() << endl;

    Point3f avgConvPoint(0, 0, 0);

    vector<Point3f> convPoints;


    for (size_t y = 0; y < final_rays.size(); y++)
    {
        for (size_t x = 0; x < final_rays[0].size(); x++)
        {
            const auto &ray = final_rays[y][x];
            Vec3b &d = dir_img.at<Vec3b>(y, x);
            if (ray.rating < 0)
            {
                d = Vec3b(0, 0, firstImage.at<Vec3b>(y, x)[0] / 4);
                continue;
            }

            // Correct ray direction: I want the final rays to be on the closest image pointing away.
            Vec3f start = ray.start + ray.dir * ((z - ray.start[2]) / ray.dir[2]);
            Vec3f dir = Vec3f(-ray.dir[0], -ray.dir[1], ray.dir[2]); // Also need to mirror the direction in order to correct the start

            // Store result with swapped axes
            result << y << ' ' << x << ' ' << start[1] << ' ' << start[0] << ' ' << dir[1] << ' ' << dir[0] << ' ' << dir[2] << endl;

            // Debug image
            d[0] = ((ray.dir[0] - min_dx) / (max_dx - min_dx)) * 255;
            d[1] = ((ray.dir[1] - min_dy) / (max_dy - min_dy)) * 255;
            d[2] = firstImage.at<Vec3b>(y, x)[0];

            // Calulate the convergence point with a other random ray
            size_t ox = rand() % final_rays[0].size();
            size_t oy = rand() % final_rays.size();

            if (abs(ox - x) > final_rays[0].size() / 4 && abs(oy - y) > final_rays.size())
            {
                const auto &other = final_rays[oy][ox];
                if (other.rating >= 0) {
                    Point3f nr = estimateConvergencePoint(ray, other);
                    /*
                    if (rand() % 10000 == 0) {
                        cout << "Px A = (" << x << ", " << y << ")" << endl;
                        cout << "Px B = (" << ox << ", " << oy << ")" << endl;
                        cout << "Ray A: s=" << ray.start << " d=" << ray.dir << endl;
                        cout << "Ray B: s=" << other.start << " d=" << other.dir << endl;

                        cout << "Ws A0=" << pixel_to_world[0][y][x] << " A1=" << pixel_to_world[1][y][x] << endl;
                        cout << "Ws B0=" << pixel_to_world[0][oy][ox] << " B1=" << pixel_to_world[1][oy][ox] << endl;
                        cout << "Conv point: " << nr << endl;

                        EstimatedRay ra;
                        ra.start = start;
                        ra.dir = dir;
                        EstimatedRay rb;
                        rb.start = other.start + other.dir * ((z - other.start[2]) / other.dir[2]);
                        rb.dir = Vec3f(-other.dir[0], -other.dir[1], other.dir[2]);

                        cout << "ra: s=" << ra.start << " d=" << ra.dir << endl;
                        cout << "rb: s=" << rb.start << " d=" << rb.dir << endl;

                        ra.start[2] = 0;
                        rb.start[2] = 0;

                        nr = estimateConvergencePoint(ra, rb);
                        cout << "Projected conv point: " << nr << endl;
                    }
                    */

                    if (isnormal(nr.x) && isnormal(nr.y) && isnormal(nr.z))
                    {
                    	avgConvPoint += nr;
                    	convPoints.push_back(nr);
                    }
                }
            }

        }
    }

    float avgConvErr = 0;
    avgConvPoint *= (1.0f / convPoints.size());

    for (const auto &c : convPoints)
    {
        float d = norm(avgConvPoint - c);
        avgConvErr += d;
    }

    avgConvErr /= convPoints.size();
    errorstr << "Convergence point is at: " << avgConvPoint << endl;
    errorstr << "Conv errors: avg=" << avgConvErr << endl;
    if (convPoints.size() > 10000)
    {
        float err;

        errorstr << "Conv point without worst 10 points: " << calculateAveragePoint(convPoints, 10, err) << " err: " << err << endl;
        errorstr << "Conv point without worst 20 points: " << calculateAveragePoint(convPoints, 10, err) << " err: " << err << endl;
        errorstr << "Conv point without worst 50 points: " << calculateAveragePoint(convPoints, 30, err) << " err: " << err << endl;
        errorstr << "Conv point without worst 100 points: " << calculateAveragePoint(convPoints, 50, err) << " err: " << err << endl;
        errorstr << "Conv point without worst 250 points: " << calculateAveragePoint(convPoints, 150, err) << " err: " << err << endl;
        errorstr << "Conv point without worst 1000 points: " << calculateAveragePoint(convPoints, 750, err) << " err: " << err << endl;
        errorstr << "Conv point without worst 10000 points: " << calculateAveragePoint(convPoints, 9000, err) << " err: " << err << endl;
    }

    imwrite(output + "_dir.png", dir_img);
    //imshow("out", dir_img);
    //while(waitKey(0) != 32);

    result.close();

    ofstream errfile;
    errfile.open (output + "_accuracy.txt", ios::out | ios::trunc);
    errfile << errorstr.str();
    cout << errorstr.str();

    return 0;
}
