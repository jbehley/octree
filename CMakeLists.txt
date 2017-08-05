project(octree)
cmake_minimum_required(VERSION 2.8)

# Release by default
# Turn on Debug with "-DCMAKE_BUILD_TYPE=Debug"
if( NOT CMAKE_BUILD_TYPE )
   SET( CMAKE_BUILD_TYPE Release )
endif()

find_package(Boost)

if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
  ADD_EXECUTABLE(example1 examples/example1.cpp)
  ADD_EXECUTABLE(example2 examples/example2.cpp)
  ADD_EXECUTABLE(example3 examples/example3.cpp)
endif()

# find gtest ...
FIND_PACKAGE(GTest)

# build gtests if possible
IF(GTEST_FOUND)
  MESSAGE("Building tests with google test framework")
	ENABLE_TESTING()
	INCLUDE_DIRECTORIES(${GTEST_INCLUDE_DIRS})
	ADD_EXECUTABLE(octree-test test/octree-test.cpp)
	TARGET_LINK_LIBRARIES(octree-test GTest::GTest GTest::Main)
	ADD_TEST(octree-test octree-test)
ENDIF()
