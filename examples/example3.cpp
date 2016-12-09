#include <iostream>
#include <cstdlib>
#include <time.h>

#include "../Octree.hpp"
#include "utils.h"

/** Example for a templated descriptor computation
 * \author behley
 */

template <class PointT>
class SimpleDescriptor
{
 public:
  SimpleDescriptor(float r, uint32_t dim) : radius_(r), dim_(dim)
  {
  }

  void compute(const PointT& query, const std::vector<PointT>& pts, const unibn::Octree<PointT>& oct,
               std::vector<float>& descriptor)
  {
    memset(&descriptor[0], 0, dim_);

    std::vector<uint32_t> neighbors;
    std::vector<float> distances;

    // template is needed to tell the compiler that radiusNeighbors is a method.
    oct.template radiusNeighbors<unibn::MaxDistance<PointT> >(query, radius_, neighbors, distances);
    for (uint32_t i = 0; i < neighbors.size(); ++i) descriptor[distances[i] / radius_ * dim_] += 1;
  }

  uint32_t dim() const
  {
    return dim_;
  }

 protected:
  float radius_;
  uint32_t dim_;
};

class Point3f
{
 public:
  Point3f(float x, float y, float z) : x(x), y(y), z(z)
  {
  }

  float x, y, z;
};

int main(int argc, char** argv)
{
  if (argc < 2)
  {
    std::cerr << "filename of point cloud missing." << std::endl;
    return -1;
  }
  std::string filename = argv[1];

  std::vector<Point3f> points;
  readPoints<Point3f>(filename, points);
  std::cout << "Read " << points.size() << " points." << std::endl;
  if (points.size() == 0)
  {
    std::cerr << "Empty point cloud." << std::endl;
    return -1;
  }

  int64_t begin, end;

  // initializing the Octree with points from point cloud.
  unibn::Octree<Point3f> octree;
  unibn::OctreeParams params;
  octree.initialize(points);

  SimpleDescriptor<Point3f> desc(0.5, 5);
  std::vector<float> values(desc.dim());
  // performing descriptor computations for each point in point cloud
  begin = clock();
  for (uint32_t i = 0; i < points.size(); ++i)
  {
    desc.compute(points[i], points, octree, values);
  }
  end = clock();
  double search_time = ((double)(end - begin) / CLOCKS_PER_SEC);
  std::cout << "Computing simple descriptor for all points took " << search_time << " seconds." << std::endl;

  octree.clear();

  return 0;
}
