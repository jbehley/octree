#include <string>
#include <map>
#include <queue>
#include <gtest/gtest.h>
#include <boost/random.hpp>


#include "../Octree.hpp"

namespace
{

class Point3f
{
  public:
    Point3f(float x, float y, float z) :
        x(x), y(y), z(z)
    {
    }

    float x, y, z;
};

// simple bruteforce search.
template<typename PointT>
class NaiveNeighborSearch
{
  public:
    void initialize(const std::vector<PointT>& points)
    {
      data_ = &points;
    }

    template<typename Distance>
    uint32_t findNeighbor(const PointT& query, float minDistance = -1.0f)
    {
      const std::vector<PointT>& pts = *data_;
      if (pts.size() == 0) return -1;

      float maxDistance = std::numeric_limits<float>::infinity();
      float sqrMinDistance = (minDistance < 0) ? minDistance : Distance::sqr(minDistance);
      uint32_t resultIndex = -1;
      for (uint32_t i = 0; i < pts.size(); ++i)
      {
        float dist = Distance::compute(query, pts[i]);
        if ((dist > sqrMinDistance) && (dist < maxDistance))
        {
          maxDistance = dist;
          resultIndex = i;
        }
      }

      return resultIndex;
    }

    template<typename Distance>
    void radiusNeighbors(const PointT& query, float radius, std::vector<uint32_t>& resultIndices)
    {
      const std::vector<PointT>& pts = *data_;
      resultIndices.clear();
      float sqrRadius = Distance::sqr(radius);

      for (uint32_t i = 0; i < pts.size(); ++i)
      {
        if (Distance::compute(query, pts[i]) < sqrRadius)
        {
          resultIndices.push_back(i);
        }
      }
    }

  protected:
    const std::vector<PointT>* data_;
};

// The fixture for testing class Foo.
class OctreeTest: public ::testing::Test
{
  public:
    typedef unibn::Octree<Point3f>::Octant Octant;

  protected:
    // helper methods to access the protected parts of octree for consistency checks.
    template<typename PointT>
    const typename unibn::Octree<PointT>::Octant* getRoot(const unibn::Octree<PointT>& oct)
    {
      return oct.root_;
    }

    template<typename PointT>
    const std::vector<uint32_t>& getSuccessors(const unibn::Octree<PointT>& oct)
    {
      return oct.successors_;
    }
};

template<typename PointT>
void randomPoints(std::vector<PointT>& pts, uint32_t N, uint32_t seed = 0)
{
  boost::mt11213b mtwister(seed);
  boost::uniform_01<> gen;
  pts.clear();
  pts.reserve(N);
  // generate N random points in [-5.0,5.0] x [-5.0,5.0] x [-5.0,5.0]...
  for (uint32_t i = 0; i < N; ++i)
  {
    float x = 10.0f * gen(mtwister) - 5.0f;
    float y = 10.0f * gen(mtwister) - 5.0f;
    float z = 10.0f * gen(mtwister) - 5.0f;

    pts.push_back(Point3f(x, y, z));
  }
}

TEST_F(OctreeTest, Initialize)
{


  uint32_t N = 1000;
  unibn::OctreeParams params;
  params.bucketSize = 16;

  unibn::Octree<Point3f> oct;

  const Octant* root = getRoot(oct);
  const std::vector<uint32_t>& successors = getSuccessors(oct);

  ASSERT_EQ(0, root);

  std::vector<Point3f> points;
  randomPoints(points, N, 1337);

  oct.initialize(points, params);

  root = getRoot(oct);

  // check first some pre-requisits.
  ASSERT_EQ(true, (root != 0));
  ASSERT_EQ(N, successors.size());

  std::vector<uint32_t> elementCount(N, 0);
  uint32_t idx = root->start;
  for (uint32_t i = 0; i < N; ++i)
  {
    ASSERT_LT(idx, N);
    ASSERT_LE(successors[idx], N);
    elementCount[idx] += 1;
    ASSERT_EQ(1, elementCount[idx]);
    idx = successors[idx];
  }

  // check that each index was found.
  for (uint32_t i = 0; i < N; ++i)
  {
    ASSERT_EQ(1, elementCount[i]);
  }

  // test if each Octant contains only points inside the octant and child octants have only real subsets of parents!
  std::queue<const Octant*> queue;
  queue.push(root);
  std::vector<int32_t> assignment(N, -1);

  while (!queue.empty())
  {
    const Octant* octant = queue.front();
    queue.pop();

    // check points.
    ASSERT_LT(octant->start, N);

    // test if each point assigned to a octant really is inside the octant.

    uint32_t idx = octant->start;
    uint32_t lastIdx = octant->start;
    for (uint32_t i = 0; i < octant->size; ++i)
    {
      float x = points[idx].x - octant->x;
      float y = points[idx].y - octant->y;
      float z = points[idx].z - octant->z;

      ASSERT_LE(std::abs(x), octant->extent);
      ASSERT_LE(std::abs(y), octant->extent);
      ASSERT_LE(std::abs(z), octant->extent);
      assignment[idx] = -1; // reset of child assignments.
      lastIdx = idx;
      idx = successors[idx];
    }
    ASSERT_EQ(octant->end, lastIdx);

    bool shouldBeLeaf = true;
    Octant* firstchild = 0;
    Octant* lastchild = 0;
    uint32_t pointSum = 0;

    for (uint32_t c = 0; c < 8; ++c)
    {
      Octant* child = octant->child[c];
      if (child == 0) continue;
      shouldBeLeaf = false;

      // child nodes should have start end intervals, which are true subsets of the parent.
      if (firstchild == 0) firstchild = child;
      // the child nodes should have intervals, where succ(e_{c-1}) == s_{c}, and \sum_c size(c) = parent size!
      if (lastchild != 0) ASSERT_EQ(child->start, successors[lastchild->end]);

      pointSum += child->size;
      lastchild = child;
      uint32_t idx = child->start;
      for (uint32_t i = 0; i < child->size; ++i)
      {
        // check if points are uniquely assigned to single child octant.
        ASSERT_EQ(-1, assignment[idx]);
        assignment[idx] = c;
        idx = successors[idx];
      }

      queue.push(child);
    }

    // consistent start/end of octant and its first and last children.
    if (firstchild != 0) ASSERT_EQ(octant->start, firstchild->start);
    if (lastchild != 0) ASSERT_EQ(octant->end, lastchild->end);

    // check leafs flag.
    ASSERT_EQ(shouldBeLeaf, octant->isLeaf);
    ASSERT_EQ((octant->size <= params.bucketSize), octant->isLeaf);

    // test if every point is assigned to a child octant.
    if (!octant->isLeaf)
    {
      ASSERT_EQ(octant->size, pointSum);
      uint32_t idx = octant->start;
      for (uint32_t i = 0; i < octant->size; ++i)
      {
        ASSERT_GT(assignment[idx], -1);
        idx = successors[idx];
      }
    }
  }
}

//TEST_F(OctreeTest, FindNeighbor)
//{
//  // compare with bruteforce search.
//  uint32_t N = 1000;
//
//  boost::mt11213b mtwister(1234);
//  boost::uniform_int<> uni_dist(0, N - 1);
//
//  std::vector<Point3f> points;
//  randomPoints(points, N, 1234);
//
//  NaiveNeighborSearch<Point3f> bruteforce;
//  bruteforce.initialize(points);
//  unibn::Octree<Point3f> octree;
//  octree.initialize(points);
//
//  for (uint32_t i = 0; i < 10; ++i)
//  {
//    uint32_t index = uni_dist(mtwister);
//    const Point3f& query = points[index];
//    // allow self-match
//    ASSERT_EQ(index, bruteforce.findNeighbor<unibn::L2Distance<Point3f> >(query));
//    ASSERT_EQ(bruteforce.findNeighbor<unibn::L2Distance<Point3f> >(query),
//        octree.findNeighbor<unibn::L2Distance<Point3f> >(query));
//
//    // disallow self-match
//    uint32_t bfneighbor = bruteforce.findNeighbor<unibn::L2Distance<Point3f> >(query, 0.0f);
//    uint32_t octneighbor = octree.findNeighbor<unibn::L2Distance<Point3f> >(query, 0.0f);
//
////    std::cout << "bf dist = " << unibn::L2Distance<Point3f>::compute(query, points[bfneighbor]) << std::endl;
////    std::cout << "oct dist = " << unibn::L2Distance<Point3f>::compute(query, points[octneighbor]) << std::endl;
//
//    ASSERT_EQ(bfneighbor, octneighbor);
//    ASSERT_EQ(bfneighbor, octree.findNeighborNoPruning<unibn::L2Distance<Point3f> >(query, 0.0f));
//
//    ASSERT_EQ(bruteforce.findNeighbor<unibn::L2Distance<Point3f> >(query, 0.3f),
//        octree.findNeighbor<unibn::L2Distance<Point3f> >(query,0.3f));
//  }
//
//}

template<typename T>
bool similarVectors(std::vector<T>& vec1, std::vector<T>& vec2)
{
  if (vec1.size() != vec2.size())
  {
    std::cout << "expected size = " << vec1.size() << ", but got size = " << vec2.size() << std::endl;
    return false;
  }

  for (uint32_t i = 0; i < vec1.size(); ++i)
  {
    bool found = false;
    for (uint32_t j = 0; j < vec2.size(); ++j)
    {
      if (vec1[i] == vec2[j])
      {
        found = true;
        break;
      }
    }
    if (!found)
    {
      std::cout << i << "-th element (" << vec1[i] << ") not found." << std::endl;
      return false;
    }
  }

  return true;
}

TEST_F(OctreeTest, RadiusNeighbors)
{
  uint32_t N = 1000;

  boost::mt11213b mtwister(1234);
  boost::uniform_int<> uni_dist(0, N - 1);

  std::vector<Point3f> points;
  randomPoints(points, N, 1234);

  NaiveNeighborSearch<Point3f> bruteforce;
  bruteforce.initialize(points);
  unibn::Octree<Point3f> octree;
  octree.initialize(points);

  float radii[4] =
  { 0.5, 1.0, 2.0, 5.0 };

  for (uint32_t r = 0; r < 4; ++r)
  {
    for (uint32_t i = 0; i < 10; ++i)
    {
      std::vector<uint32_t> neighborsBruteforce;
      std::vector<uint32_t> neighborsOctree;

      const Point3f& query = points[uni_dist(mtwister)];

      bruteforce.radiusNeighbors<unibn::L2Distance<Point3f> >(query, radii[r], neighborsBruteforce);
      octree.radiusNeighbors<unibn::L2Distance<Point3f> >(query, radii[r], neighborsOctree);
      ASSERT_EQ(true, similarVectors(neighborsBruteforce, neighborsOctree));

      bruteforce.radiusNeighbors<unibn::L1Distance<Point3f> >(query, radii[r], neighborsBruteforce);
      octree.radiusNeighbors<unibn::L1Distance<Point3f> >(query, radii[r], neighborsOctree);

      ASSERT_EQ(true, similarVectors(neighborsBruteforce, neighborsOctree));

      bruteforce.radiusNeighbors<unibn::MaxDistance<Point3f> >(query, radii[r], neighborsBruteforce);
      octree.radiusNeighbors<unibn::MaxDistance<Point3f> >(query, radii[r], neighborsOctree);

      ASSERT_EQ(true, similarVectors(neighborsBruteforce, neighborsOctree));
    }
  }
}

}

int main(int argc, char **argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
