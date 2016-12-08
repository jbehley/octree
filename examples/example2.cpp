#include <iostream>
#include <cstdlib>
#include <time.h>

#include "../Octree.hpp"
#include "utils.h"

/** Example 2: Same as first example, but with different point data structure, showing the
 *  variability of the Octree with access traits.
 *
 * \author behley
 */

class CustomPoint
{
 public:
  CustomPoint(float x, float y, float z) : x(x), y(y), z(z)
  {
  }

  float getX() const
  {
    return x;
  }
  float getY() const
  {
    return y;
  }
  float getZ() const
  {
    return z;
  }

 protected:
  float x, y, z;
};

// access traits for CustomPoint
namespace unibn
{
namespace traits
{

template <>
struct access<CustomPoint, 0>
{
  static float get(const CustomPoint& p)
  {
    return p.getX();
  }
};

template <>
struct access<CustomPoint, 1>
{
  static float get(const CustomPoint& p)
  {
    return p.getY();
  }
};

template <>
struct access<CustomPoint, 2>
{
  static float get(const CustomPoint& p)
  {
    return p.getZ();
  }
};
}
}

int main(int argc, char** argv)
{
  if (argc < 2)
  {
    std::cerr << "filename of point cloud missing." << std::endl;
    return -1;
  }
  std::string filename = argv[1];

  std::vector<CustomPoint> points;
  readPoints<CustomPoint>(filename, points);
  std::cout << "Read " << points.size() << " points." << std::endl;
  if (points.size() == 0)
  {
    std::cerr << "Empty point cloud." << std::endl;
    return -1;
  }

  int64_t begin, end;

  // initializing the Octree with points from point cloud.
  unibn::Octree<CustomPoint> octree;
  unibn::OctreeParams params;
  octree.initialize(points);

  // radiusNeighbors returns indexes to neighboring points.
  std::vector<uint32_t> results;
  const CustomPoint& q = points[0];
  octree.radiusNeighbors<unibn::L2Distance<CustomPoint> >(q, 0.2f, results);
  std::cout << results.size() << " radius neighbors (r = 0.2m) found for (" << q.getX() << ", " << q.getY() << ","
            << q.getZ() << ")" << std::endl;
  for (uint32_t i = 0; i < results.size(); ++i)
  {
    const CustomPoint& p = points[results[i]];
    std::cout << "  " << results[i] << ": (" << p.getX() << ", " << p.getY() << ", " << p.getZ() << ") => "
              << std::sqrt(unibn::L2Distance<CustomPoint>::compute(p, q)) << std::endl;
  }

  // performing queries for each point in point cloud
  begin = clock();
  for (uint32_t i = 0; i < points.size(); ++i)
  {
    octree.radiusNeighbors<unibn::L2Distance<CustomPoint> >(points[i], 0.5f, results);
  }
  end = clock();
  double search_time = ((double)(end - begin) / CLOCKS_PER_SEC);
  std::cout << "Searching for all radius neighbors (r = 0.5m) took " << search_time << " seconds." << std::endl;

  octree.clear();

  return 0;
}
