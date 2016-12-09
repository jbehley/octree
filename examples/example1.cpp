#include <iostream>
#include <cstdlib>
#include <time.h>

#include "../Octree.hpp"
#include "utils.h"

/** Example 1: Searching radius neighbors with default access by public x,y,z variables.
 *
 * \author behley
 */

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

  // radiusNeighbors returns indexes to neighboring points.
  std::vector<uint32_t> results;
  const Point3f& q = points[0];
  octree.radiusNeighbors<unibn::L2Distance<Point3f> >(q, 0.2f, results);
  std::cout << results.size() << " radius neighbors (r = 0.2m) found for (" << q.x << ", " << q.y << "," << q.z << ")"
            << std::endl;
  for (uint32_t i = 0; i < results.size(); ++i)
  {
    const Point3f& p = points[results[i]];
    std::cout << "  " << results[i] << ": (" << p.x << ", " << p.y << ", " << p.z << ") => "
              << std::sqrt(unibn::L2Distance<Point3f>::compute(p, q)) << std::endl;
  }

  // performing queries for each point in point cloud
  begin = clock();
  for (uint32_t i = 0; i < points.size(); ++i)
  {
    octree.radiusNeighbors<unibn::L2Distance<Point3f> >(points[i], 0.5f, results);
  }
  end = clock();
  double search_time = ((double)(end - begin) / CLOCKS_PER_SEC);
  std::cout << "Searching for all radius neighbors (r = 0.5m) took " << search_time << " seconds." << std::endl;

  octree.clear();

  return 0;
}
