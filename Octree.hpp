#ifndef OCTREE_H_
#define OCTREE_H_

// Copyright (c) 2015 Jens Behley, University of Bonn
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to
// deal in the Software without restriction, including without limitation the
// rights  to use, copy, modify, merge, publish, distribute, sublicense, and/or
// sell copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
// IN THE SOFTWARE.

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>   // memset.
#include <iostream>  // NOLINT
#include <limits>
#include <stdint.h>
#include <vector>

// needed for gtest access to protected/private members ...
/*namespace {
class OctreeTest;
}*/

namespace unibn {

/**
 * Some traits to access coordinates regardless of the specific implementation
 *of point
 * inspired by boost.geometry, which needs to be implemented by new points.
 *
 */
namespace traits {

template <typename PointT, int D>
struct access {};

template <class PointT>
struct access<PointT, 0> {
  static float get(const PointT& p) { return p.x; }
};

template <class PointT>
struct access<PointT, 1> {
  static float get(const PointT& p) { return p.y; }
};

template <class PointT>
struct access<PointT, 2> {
  static float get(const PointT& p) { return p.z; }
};
}  // namespace traits

/** convenience function for access of point coordinates **/
template <int D, typename PointT>
inline float get(const PointT& p) {
  return traits::access<PointT, D>::get(p);
}

/** convenience functions added to this library to
 * do some useful operation on Point **/
template <typename PointT>
PointT cross(PointT u, PointT v) {
    PointT cross;
    cross.x = u.y*v.z - u.z*v.y;
    cross.y = u.z*v.x - u.x*v.z;
    cross.z = u.x*v.y - u.y*v.x;
    return cross;
}
template <typename PointT>
PointT sub(PointT u, PointT v) {
    PointT sub;
    sub.x = u.x - v.x;
    sub.y = u.y - v.y;
    sub.z = u.z - v.z;
    return sub;
}
template <typename PointT>
float dot(PointT u, PointT v) {
    float dot;
    dot = u.x * v.x + u.y * v.y + u.z * v.z;
    return dot;
}

/**
 * Some generic distances: Manhattan, (squared) Euclidean, and Maximum distance.
 *
 * A Distance has to implement the methods
 * 1. compute of two points p and q to compute and return the distance between
 *two points, and
 * 2. norm of x,y,z coordinates to compute and return the norm of a point p =
 *(x,y,z)
 * 3. sqr and sqrt of value to compute the correct radius if a comparison is
 *performed using squared norms (see L2Distance)...
 */
template <typename PointT>
struct L1Distance {
  static inline float compute(const PointT& p, const PointT& q) {
    float diff1 = get<0>(p) - get<0>(q);
    float diff2 = get<1>(p) - get<1>(q);
    float diff3 = get<2>(p) - get<2>(q);

    return std::abs(diff1) + std::abs(diff2) + std::abs(diff3);
  }

  static inline float norm(float x, float y, float z) {
    return std::abs(x) + std::abs(y) + std::abs(z);
  }

  static inline float sqr(float r) { return r; }

  static inline float sqrt(float r) { return r; }
};

template <typename PointT>
struct L2Distance {
  static inline float compute(const PointT& p, const PointT& q) {
    float diff1 = get<0>(p) - get<0>(q);
    float diff2 = get<1>(p) - get<1>(q);
    float diff3 = get<2>(p) - get<2>(q);

    return std::pow(diff1, 2) + std::pow(diff2, 2) + std::pow(diff3, 2);
  }

  static inline float norm(float x, float y, float z) {
    return std::pow(x, 2) + std::pow(y, 2) + std::pow(z, 2);
  }

  static inline float sqr(float r) { return r * r; }

  static inline float sqrt(float r) { return std::sqrt(r); }
};

template <typename PointT>
struct MaxDistance {
  static inline float compute(const PointT& p, const PointT& q) {
    float diff1 = std::abs(get<0>(p) - get<0>(q));
    float diff2 = std::abs(get<1>(p) - get<1>(q));
    float diff3 = std::abs(get<2>(p) - get<2>(q));

    float maximum = diff1;
    if (diff2 > maximum) maximum = diff2;
    if (diff3 > maximum) maximum = diff3;

    return maximum;
  }

  static inline float norm(float x, float y, float z) {
    float maximum = x;
    if (y > maximum) maximum = y;
    if (z > maximum) maximum = z;
    return maximum;
  }

  static inline float sqr(float r) { return r; }

  static inline float sqrt(float r) { return r; }
};

struct OctreeParams {
 public:
  explicit OctreeParams(uint32_t bucketSize = 32, bool copyPoints = false)
      : bucketSize(bucketSize), copyPoints(copyPoints) {}
  uint32_t bucketSize;
  bool copyPoints;
};

/** \brief Index-based Octree implementation offering different queries and
 *insertion/removal of points.
 *
 * The index-based Octree uses a successor relation and a startIndex in each
 *Octant to improve runtime
 * performance for radius queries. The efficient storage of the points by
 *relinking list elements
 * bases on the insight that children of an Octant contain disjoint subsets of
 *points inside the Octant and
 * that we can reorganize the points such that we get an continuous single
 *connect list that we can use to
 * store in each octant the start of this list.
 *
 * Special about the implementation is that it allows to search for neighbors
 *with arbitrary p-norms, which
 * distinguishes it from most other Octree implementations.
 *
 * We decided to implement the Octree using a template for points and
 *containers. The container must have an
 * operator[], which allows to access the points, and a size() member function,
 *which allows to get the size of the
 * container. For the points, we used an access trait to access the coordinates
 *inspired by boost.geometry.
 * The implementation already provides a general access trait, which expects to
 *have public member variables x,y,z.
 *
 * f you use the implementation or ideas from the corresponding paper in your
 *academic work, it would be nice if you
 * cite the corresponding paper:
 *
 *    J. Behley, V. Steinhage, A.B. Cremers. Efficient Radius Neighbor Search in
 *Three-dimensional Point Clouds,
 *    Proc. of the IEEE International Conference on Robotics and Automation
 *(ICRA), 2015.
 *
 * In future, we might add also other neighbor queries and implement the removal
 *and adding of points.
 *
 * \version 0.1-icra
 *
 * \author behley
 */

template <typename PointT, typename ContainerT = std::vector<PointT> >
class Octree {
 public:
  Octree();
  ~Octree();

  /** \brief initialize octree with all points **/
  void initialize(const ContainerT& pts,
                  const OctreeParams& params = OctreeParams());

  /** \brief initialize octree only from pts that are inside indexes. **/
  void initialize(const ContainerT& pts, const std::vector<uint32_t>& indexes,
                  const OctreeParams& params = OctreeParams());

  /** \brief remove all data inside the octree. **/
  void clear();

  /** \brief function added to this library to expand the octree to include new_pt. **/
  template <typename Distance>
  void expand(const PointT& new_pt, const ContainerT& pts,
              const OctreeParams& params = OctreeParams());

  /** \brief radius neighbor queries where radius determines the maximal radius
   * of reported indices of points in resultIndices **/
  template <typename Distance>
  void radiusNeighbors(const PointT& query, float radius,
                       std::vector<uint32_t>* const resultIndices) const;

  /** \brief radius neighbor queries with explicit (squared) distance
   * computation. **/
  template <typename Distance>
  void radiusNeighbors(const PointT& query, float radius,
                       std::vector<uint32_t>* const resultIndices,
                       std::vector<float>* const distances) const;

  /** \brief function added to this library to perform
   * nearest neighbor queries. The index of the point found is in resultIndice **/
  template <typename Distance>
  void nearestNeighbor(const PointT& query, uint32_t* const resultIndice) const;

 protected:
  class Octant {
   public:
    Octant();
    ~Octant();

    bool isLeaf;

    // bounding box of the octant needed for overlap and contains tests...
    float x, y, z;  // center
    float extent;   // half of side-length

    uint32_t start, end;  // start and end in succ_
    uint32_t size;        // number of points

    Octant* child[8];
  };

  // not copyable, not assignable ...
  Octree(Octree&);
  Octree& operator=(const Octree& oct);

  /**
   * \brief creation of an octant using the elements starting at startIdx.
   *
   * The method reorders the index such that all points are correctly linked to
   *successors belonging
   * to the same octant.
   *
   * \param x,y,z           center coordinates of octant
   * \param extent          extent of octant
   * \param startIdx        first index of points inside octant
   * \param endIdx          last index of points inside octant
   * \param size            number of points in octant
   *
   * \return  octant with children nodes.
   */
  Octant* createOctant(float x, float y, float z, float extent,
                       uint32_t startIdx, uint32_t endIdx, uint32_t size);

  /** \brief function added to this library to test if point is inside octant cell.
   * To do so the 8 corners of the octant cell are used.
   * These corners are defined as follow :
   *
   *         7         6
   *         *--------*
   *        /|       /|
   *       / |      / |
   *     3*--------*2 |
   *      | 4*-----|--*5
   *      | /      | /
   *      |/       |/
   *      *--------*
   *      0        1
   *
   * @param point         point to test
   * @param octant_x      x coordinate of center of the octant
   * @param octant_y      y coordinate of center of the octant
   * @param octant_z      z coordinate of center of the octant
   * @param octant_extent extent of the octant
   *
   * @return true if point is inside octant cell, false otherwise.
   */
  template <typename Distance>
  bool insideOctantCell(const PointT& point, float octant_x, float octant_y,
                        float octant_z, float octant_extent);

  template <typename Distance>
  void radiusNeighbors(const Octant* octant, const PointT& query, float radius,
                       float sqrRadius,
                       std::vector<uint32_t>* const resultIndices) const;

  template <typename Distance>
  void radiusNeighbors(const Octant* octant, const PointT& query, float radius,
                       float sqrRadius,
                       std::vector<uint32_t>* const resultIndices,
                       std::vector<float>* const distances) const;

  template <typename Distance>
  void nearestNeighbor(const Octant* octant, const PointT& query,
                       float* const radius, float* const closest_distance,
                       uint32_t* const resultIndice) const;

  /** \brief test if search ball S(q,r) overlaps with octant
   *
   * @param query   query point
   * @param radius  "squared" radius
   * @param o       pointer to octant
   *
   * @return true, if search ball overlaps with octant, false otherwise.
   */
  template <typename Distance>
  static bool overlaps(const PointT& query, float radius, float sqRadius,
                       const Octant* o);

  /** \brief test if search ball S(q,r) contains octant
   *
   * @param query    query point
   * @param sqRadius "squared" radius
   * @param octant   pointer to octant
   *
   * @return true, if search ball overlaps with octant, false otherwise.
   */
  template <typename Distance>
  static bool contains(const PointT& query, float sqRadius,
                       const Octant* octant);

  OctreeParams params_;
  Octant* root_;
  const ContainerT* data_;

  std::vector<uint32_t> successors_;  // single connected list of next point
                                      // indices...

  // friend class ::OctreeTest;
};

template <typename PointT, typename ContainerT>
Octree<PointT, ContainerT>::Octant::Octant()
    : isLeaf(true),
      x(0.0f),
      y(0.0f),
      z(0.0f),
      extent(0.0f),
      start(0),
      end(0),
      size(0) {
  memset(&child, 0, 8 * sizeof(Octant*));
}

template <typename PointT, typename ContainerT>
Octree<PointT, ContainerT>::Octant::~Octant() {
  for (uint32_t i = 0; i < 8; ++i) delete child[i];
}

template <typename PointT, typename ContainerT>
Octree<PointT, ContainerT>::Octree()
    : root_(0), data_(0) {}

template <typename PointT, typename ContainerT>
Octree<PointT, ContainerT>::~Octree() {
  delete root_;
  if (params_.copyPoints) delete data_;
}

template <typename PointT, typename ContainerT>
void Octree<PointT, ContainerT>::initialize(const ContainerT& pts,
                                            const OctreeParams& params) {
  clear();
  params_ = params;

  if (params_.copyPoints)
    data_ = new ContainerT(pts);
  else
    data_ = &pts;

  const uint32_t N = pts.size();
  successors_ = std::vector<uint32_t>(N);

  // determine axis-aligned bounding box.
  float min[3], max[3];
  min[0] = get<0>(pts[0]);
  min[1] = get<1>(pts[0]);
  min[2] = get<2>(pts[0]);
  max[0] = min[0];
  max[1] = min[1];
  max[2] = min[2];

  for (uint32_t i = 0; i < N; ++i) {
    // initially each element links simply to the following element.
    successors_[i] = i + 1;

    const PointT& p = pts[i];

    if (get<0>(p) < min[0]) min[0] = get<0>(p);
    if (get<1>(p) < min[1]) min[1] = get<1>(p);
    if (get<2>(p) < min[2]) min[2] = get<2>(p);
    if (get<0>(p) > max[0]) max[0] = get<0>(p);
    if (get<1>(p) > max[1]) max[1] = get<1>(p);
    if (get<2>(p) > max[2]) max[2] = get<2>(p);
  }

  float ctr[3] = {min[0], min[1], min[2]};

  float maxextent = 0.5f * (max[0] - min[0]);
  ctr[0] += maxextent;
  for (uint32_t i = 1; i < 3; ++i) {
    float extent = 0.5f * (max[i] - min[i]);
    ctr[i] += extent;
    if (extent > maxextent) maxextent = extent;
  }

  root_ = createOctant(ctr[0], ctr[1], ctr[2], maxextent, 0, N - 1, N);
}

template <typename PointT, typename ContainerT>
void Octree<PointT, ContainerT>::initialize(
    const ContainerT& pts, const std::vector<uint32_t>& indexes,
    const OctreeParams& params) {
  clear();
  params_ = params;
  const uint32_t N = pts.size();
  successors_ = std::vector<uint32_t>(N);

  if (indexes.size() == 0) return;

  // determine axis-aligned bounding box.
  uint32_t lastIdx = indexes[0];
  float min[3], max[3];
  min[0] = get<0>(pts[lastIdx]);
  min[1] = get<1>(pts[lastIdx]);
  min[2] = get<2>(pts[lastIdx]);
  max[0] = min[0];
  max[1] = min[1];
  max[2] = min[2];

  for (uint32_t i = 1; i < indexes.size(); ++i) {
    uint32_t idx = indexes[i];
    // initially each element links simply to the following element.
    successors_[lastIdx] = idx;

    const PointT& p = pts[idx];

    if (get<0>(p) < min[0]) min[0] = get<0>(p);
    if (get<1>(p) < min[1]) min[1] = get<1>(p);
    if (get<2>(p) < min[2]) min[2] = get<2>(p);
    if (get<0>(p) > max[0]) max[0] = get<0>(p);
    if (get<1>(p) > max[1]) max[1] = get<1>(p);
    if (get<2>(p) > max[2]) max[2] = get<2>(p);

    lastIdx = idx;
  }

  float ctr[3] = {min[0], min[1], min[2]};

  float maxextent = 0.5f * (max[0] - min[0]);
  ctr[0] += maxextent;
  for (uint32_t i = 1; i < 3; ++i) {
    float extent = 0.5f * (max[i] - min[i]);
    ctr[i] += extent;
    if (extent > maxextent) maxextent = extent;
  }

  root_ = createOctant(ctr[0], ctr[1], ctr[2], maxextent, indexes[0], lastIdx,
                       indexes.size());
}

template <typename PointT, typename ContainerT>
void Octree<PointT, ContainerT>::clear() {
  delete root_;
  if (params_.copyPoints) delete data_;
  root_ = 0;
  data_ = 0;
  successors_.clear();
}

template <typename PointT, typename ContainerT>
template<typename Distance>
void Octree<PointT, ContainerT>::expand(const PointT& new_pt, const ContainerT& pts,
                                        const OctreeParams& params) {
  float prev_root_x = root_->x;
  float prev_root_y = root_->y;
  float prev_root_z = root_->z;
  float prev_root_extent = root_->extent;
  clear();
  params_ = params;
  if (params_.copyPoints)
    data_ = new ContainerT(pts);
  else
    data_ = &pts;
  const uint32_t N = pts.size();
  successors_ = std::vector<uint32_t>(N);
  for (uint32_t i = 0; i < N; ++i)
    successors_[i] = i + 1;
  // expand the root cell until new_point is located inside
  while (!insideOctantCell<Distance>(new_pt, prev_root_x, prev_root_y, prev_root_z, prev_root_extent)) {
    // chose center of new root as closest corner of current root to new_point
    std::vector<PointT> corners;
    std::vector<float> distances_to_corners;
    corners.push_back(PointT(prev_root_x - prev_root_extent, prev_root_y - prev_root_extent, prev_root_z - prev_root_extent));
    corners.push_back(PointT(prev_root_x + prev_root_extent, prev_root_y - prev_root_extent, prev_root_z - prev_root_extent));
    corners.push_back(PointT(prev_root_x + prev_root_extent, prev_root_y + prev_root_extent, prev_root_z - prev_root_extent));
    corners.push_back(PointT(prev_root_x - prev_root_extent, prev_root_y + prev_root_extent, prev_root_z - prev_root_extent));
    corners.push_back(PointT(prev_root_x - prev_root_extent, prev_root_y - prev_root_extent, prev_root_z + prev_root_extent));
    corners.push_back(PointT(prev_root_x + prev_root_extent, prev_root_y - prev_root_extent, prev_root_z + prev_root_extent));
    corners.push_back(PointT(prev_root_x + prev_root_extent, prev_root_y + prev_root_extent, prev_root_z + prev_root_extent));
    corners.push_back(PointT(prev_root_x - prev_root_extent, prev_root_y + prev_root_extent, prev_root_z + prev_root_extent));
    for (unsigned int i = 0; i < corners.size(); ++i)
      distances_to_corners.push_back(Distance::compute(new_pt, corners[i]));
    std::vector<float>::iterator closest_corner_it = std::min_element(distances_to_corners.begin(), distances_to_corners.end());
    int closest_corner_index = std::distance(distances_to_corners.begin(), closest_corner_it) - 1;
    prev_root_x = corners[closest_corner_index].x;
    prev_root_y = corners[closest_corner_index].y;
    prev_root_z = corners[closest_corner_index].z;
    prev_root_extent *= 2;
  }
  root_ = createOctant(prev_root_x, prev_root_y, prev_root_z, prev_root_extent, 0, N - 1, N);
}

template <typename PointT, typename ContainerT>
template<typename Distance>
bool Octree<PointT, ContainerT>::insideOctantCell(const PointT& point, float octant_x, float octant_y, float octant_z, float octant_extent) {
  std::vector<PointT> corners;
  PointT above_normal, below_normal, left_normal, right_normal, front_normal, back_normal;
  bool is_above, is_below, is_left, is_right, is_front, is_back;
  corners.push_back(PointT(octant_x - octant_extent, octant_y - octant_extent, octant_z - octant_extent));
  corners.push_back(PointT(octant_x + octant_extent, octant_y - octant_extent, octant_z - octant_extent));
  corners.push_back(PointT(octant_x + octant_extent, octant_y + octant_extent, octant_z - octant_extent));
  corners.push_back(PointT(octant_x - octant_extent, octant_y + octant_extent, octant_z - octant_extent));
  corners.push_back(PointT(octant_x - octant_extent, octant_y - octant_extent, octant_z + octant_extent));
  corners.push_back(PointT(octant_x + octant_extent, octant_y - octant_extent, octant_z + octant_extent));
  corners.push_back(PointT(octant_x + octant_extent, octant_y + octant_extent, octant_z + octant_extent));
  corners.push_back(PointT(octant_x - octant_extent, octant_y + octant_extent, octant_z + octant_extent));

  above_normal = cross(sub(corners[3], corners[2]), sub(corners[6], corners[2]));
  below_normal = cross(sub(corners[5], corners[1]), sub(corners[0], corners[1]));
  left_normal  = cross(sub(corners[7], corners[4]), sub(corners[0], corners[4]));
  right_normal = cross(sub(corners[1], corners[5]), sub(corners[6], corners[5]));
  front_normal = cross(sub(corners[0], corners[1]), sub(corners[2], corners[1]));
  back_normal  = cross(sub(corners[6], corners[5]), sub(corners[4], corners[5]));

  is_above = dot(sub(point, corners[3]), above_normal) >= 0.;
  is_below = dot(sub(point, corners[0]), below_normal) >= 0.;
  is_left  = dot(sub(point, corners[0]), left_normal)  >= 0.;
  is_right = dot(sub(point, corners[1]), right_normal) >= 0.;
  is_front = dot(sub(point, corners[0]), front_normal) >= 0.;
  is_back  = dot(sub(point, corners[4]), back_normal)  >= 0.;
  if (!is_above && !is_below && !is_left && !is_right && !is_front && !is_back)
    return true;

  return false;
}

template <typename PointT, typename ContainerT>
typename Octree<PointT, ContainerT>::Octant*
Octree<PointT, ContainerT>::createOctant(float x, float y, float z,
                                         float extent, uint32_t startIdx,
                                         uint32_t endIdx, uint32_t size) {
  // For a leaf we don't have to change anything; points are already correctly
  // linked or correctly reordered.
  Octant* octant = new Octant;

  octant->isLeaf = true;

  octant->x = x;
  octant->y = y;
  octant->z = z;
  octant->extent = extent;

  octant->start = startIdx;
  octant->end = endIdx;
  octant->size = size;

  static const float factor[] = {-0.5f, 0.5f};
  // subdivide subset of points and re-link points according to Morton codes
  if (size > params_.bucketSize) {
    octant->isLeaf = false;

    const ContainerT& points = *data_;
    std::vector<uint32_t> childStarts(8, 0);
    std::vector<uint32_t> childEnds(8, 0);
    std::vector<uint32_t> childSizes(8, 0);

    // re-link disjoint child subsets...
    uint32_t idx = startIdx;

    for (uint32_t i = 0; i < size; ++i) {
      const PointT& p = points[idx];

      // determine Morton code for each point...
      uint32_t mortonCode = 0;
      if (get<0>(p) > x) mortonCode |= 1;
      if (get<1>(p) > y) mortonCode |= 2;
      if (get<2>(p) > z) mortonCode |= 4;

      // set child starts and update successors...
      if (childSizes[mortonCode] == 0)
        childStarts[mortonCode] = idx;
      else
        successors_[childEnds[mortonCode]] = idx;
      childSizes[mortonCode] += 1;

      childEnds[mortonCode] = idx;
      idx = successors_[idx];
    }

    // now, we can create the child nodes...
    float childExtent = 0.5f * extent;
    bool firsttime = true;
    uint32_t lastChildIdx = 0;
    for (uint32_t i = 0; i < 8; ++i) {
      if (childSizes[i] == 0) continue;

      float childX = x + factor[(i & 1) > 0] * extent;
      float childY = y + factor[(i & 2) > 0] * extent;
      float childZ = z + factor[(i & 4) > 0] * extent;

      octant->child[i] =
          createOctant(childX, childY, childZ, childExtent, childStarts[i],
                       childEnds[i], childSizes[i]);

      if (firsttime)
        octant->start = octant->child[i]->start;
      else
        successors_[octant->child[lastChildIdx]->end] =
            octant->child[i]->start;  // we have to ensure that also the child
                                      // ends link to the next child start.

      lastChildIdx = i;
      octant->end = octant->child[i]->end;
      firsttime = false;
    }
  }
  return octant;
}

template <typename PointT, typename ContainerT>
template <typename Distance>
void Octree<PointT, ContainerT>::radiusNeighbors(
    const Octant* octant, const PointT& query, float radius, float sqrRadius,
    std::vector<uint32_t>* const resultIndices) const {
  const ContainerT& points = *data_;

  // if search ball S(q,r) contains octant, simply add point indexes.
  if (contains<Distance>(query, sqrRadius, octant)) {
    uint32_t idx = octant->start;
    for (uint32_t i = 0; i < octant->size; ++i) {
      resultIndices->push_back(idx);
      idx = successors_[idx];
    }

    return;  // early pruning.
  }

  if (octant->isLeaf) {
    uint32_t idx = octant->start;
    for (uint32_t i = 0; i < octant->size; ++i) {
      const PointT& p = points[idx];
      float dist = Distance::compute(query, p);
      if (dist < sqrRadius) resultIndices->push_back(idx);
      idx = successors_[idx];
    }

    return;
  }

  // check whether child nodes are in range.
  for (uint32_t c = 0; c < 8; ++c) {
    if (octant->child[c] == 0) continue;
    if (!overlaps<Distance>(query, radius, sqrRadius, octant->child[c]))
      continue;
    radiusNeighbors<Distance>(octant->child[c], query, radius, sqrRadius,
                              resultIndices);
  }
}

template <typename PointT, typename ContainerT>
template <typename Distance>
void Octree<PointT, ContainerT>::radiusNeighbors(
    const Octant* octant, const PointT& query, float radius, float sqrRadius,
    std::vector<uint32_t>* const resultIndices,
    std::vector<float>* const distances) const {
  const ContainerT& points = *data_;

  // if search ball S(q,r) contains octant, simply add point indexes and compute
  // squared distances.
  if (contains<Distance>(query, sqrRadius, octant)) {
    uint32_t idx = octant->start;
    for (uint32_t i = 0; i < octant->size; ++i) {
      resultIndices->push_back(idx);
      distances->push_back(Distance::compute(query, points[idx]));
      idx = successors_[idx];
    }

    return;  // early pruning.
  }

  if (octant->isLeaf) {
    uint32_t idx = octant->start;
    for (uint32_t i = 0; i < octant->size; ++i) {
      const PointT& p = points[idx];
      float dist = Distance::compute(query, p);
      if (dist < sqrRadius) {
        resultIndices->push_back(idx);
        distances->push_back(dist);
      }
      idx = successors_[idx];
    }

    return;
  }

  // check whether child nodes are in range.
  for (uint32_t c = 0; c < 8; ++c) {
    if (octant->child[c] == 0) continue;
    if (!overlaps<Distance>(query, radius, sqrRadius, octant->child[c]))
      continue;
    radiusNeighbors<Distance>(octant->child[c], query, radius, sqrRadius,
                              resultIndices, distances);
  }
}

template <typename PointT, typename ContainerT>
template <typename Distance>
void Octree<PointT, ContainerT>::radiusNeighbors(
    const PointT& query, float radius,
    std::vector<uint32_t>* const resultIndices) const {
  resultIndices->clear();
  if (root_ == 0) return;

  float sqrRadius = Distance::sqr(radius);  // "squared" radius
  radiusNeighbors<Distance>(root_, query, radius, sqrRadius, resultIndices);
}

template <typename PointT, typename ContainerT>
template <typename Distance>
void Octree<PointT, ContainerT>::radiusNeighbors(
    const PointT& query, float radius,
    std::vector<uint32_t>* const resultIndices,
    std::vector<float>* const distances) const {
  resultIndices->clear();
  distances->clear();
  if (root_ == 0) return;

  float sqrRadius = Distance::sqr(radius);  // "squared" radius
  radiusNeighbors<Distance>(root_, query, radius, sqrRadius, resultIndices,
                            distances);
}

template <typename PointT, typename ContainerT>
template <typename Distance>
void Octree<PointT, ContainerT>::nearestNeighbor(
    const Octant* octant, const PointT& query, float* const radius,
    float* const closest_distance, uint32_t* const resultIndice) const {
  const ContainerT& points = *data_;

  // if search ball S(q,r) contains octant, simply pick closest point index.
  if (contains<Distance>(query, Distance::sqr(*radius), octant)) {
    uint32_t idx = octant->start;
    for (uint32_t i = 0; i < octant->size; ++i) {
      const PointT& p = points[idx];
      float dist = Distance::compute(query, p);
      if (dist < *closest_distance) {
        *resultIndice = idx;
        *closest_distance = dist;
        *radius = dist;
      }
      idx = successors_[idx];
    }
    return;  // early pruning.
  }

  if (octant->isLeaf) {
    uint32_t idx = octant->start;
    for (uint32_t i = 0; i < octant->size; ++i) {
      const PointT& p = points[idx];
      float dist = Distance::compute(query, p);
      if (dist < Distance::sqr(*radius)) {
        *resultIndice = idx;
        *radius = dist;
      }
      idx = successors_[idx];
    }
    return;
  }

  // check whether child nodes are in range.
  for (uint32_t c = 0; c < 8; ++c) {
    if (octant->child[c] == 0) continue;
    if (!overlaps<Distance>(query, *radius, Distance::sqr(*radius),
                            octant->child[c]))
      continue;
    nearestNeighbor<Distance>(octant->child[c], query, radius, closest_distance,
                              resultIndice);
  }
}

template <typename PointT, typename ContainerT>
template <typename Distance>
void Octree<PointT, ContainerT>::nearestNeighbor(
    const PointT& query, uint32_t* const resultIndice) const {
  *resultIndice = 0;
  if (root_ == 0) return;

  float radius = std::numeric_limits<double>::max();
  float closest_distance = std::numeric_limits<double>::max();
  nearestNeighbor<Distance>(root_, query, &radius, &closest_distance,
                            resultIndice);
}

template <typename PointT, typename ContainerT>
template <typename Distance>
bool Octree<PointT, ContainerT>::overlaps(const PointT& query, float radius,
                                          float sqRadius, const Octant* o) {
  // we exploit the symmetry to reduce the test to testing if its inside the
  // Minkowski sum around the positive quadrant.
  float x = get<0>(query) - o->x;
  float y = get<1>(query) - o->y;
  float z = get<2>(query) - o->z;

  x = std::abs(x);
  y = std::abs(y);
  z = std::abs(z);

  // (1) checking the line region.
  float maxdist = radius + o->extent;

  // a. completely outside, since q' is outside the relevant area.
  if (x > maxdist || y > maxdist || z > maxdist) return false;

  // b. inside the line region, one of the coordinates is inside the square.
  if (x < o->extent || y < o->extent || z < o->extent) return true;

  // (2) checking the corner region...
  x -= o->extent;
  y -= o->extent;
  z -= o->extent;

  return (Distance::norm(x, y, z) < sqRadius);
}

template <typename PointT, typename ContainerT>
template <typename Distance>
bool Octree<PointT, ContainerT>::contains(const PointT& query, float sqRadius,
                                          const Octant* o) {
  // we exploit the symmetry to reduce the test to test
  // whether the farthest corner is inside the search ball.
  float x = get<0>(query) - o->x;
  float y = get<1>(query) - o->y;
  float z = get<2>(query) - o->z;

  x = std::abs(x);
  y = std::abs(y);
  z = std::abs(z);
  // reminder: (x, y, z) - (-e, -e, -e) = (x, y, z) + (e, e, e)
  x += o->extent;
  y += o->extent;
  z += o->extent;

  return (Distance::norm(x, y, z) < sqRadius);
}
}  // namespace unibn

#endif  // OCTREE_H_

