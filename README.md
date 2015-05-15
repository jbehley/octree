
# Efficient Radius Neighbor Search in Three-dimensional Point Clouds

This repository provides the Octree implementation of the paper "Efficient Radius Neighbor Search in Three-dimensional Point Clouds" by J. Behley, V. Steinhage, and A.B. Cremers, University of Bonn presented at the IEEE International Conference on Robotics and Automation (ICRA), 2015.


## Features
- Fast radius neighbor search for three-dimensional point clouds.
- Fully templated for maximal flexibility to support arbitrary point representations & containers
- Supports arbitrary p-norms: L1, L2 and Maximum norm included.

## Building the examples & tests

The octree itself has no dependencies. 
However, for compiling the examples, you need [CMake](http://www.cmake.org/) and [Boost C++ library](http://www.boost.org/).
For building the examples you have to first build the project:

```bash
mkdir build
cd build
cmake ..
make
```

To run the examples, you need some point cloud data:

```bash
wget http://www.iai.uni-bonn.de/~behley/data/wachtberg_folds.zip
unzip wachtberg_folds.zip -d data
```

Now you can run the examples:

```bash
./example1 data/scan_001_points.dat
```

which perform some queries and demonstrate the flexibility of our Octree implementation to handle different implementations of points.

The different examples show some use cases of the octree. `example1` demonstrates the general usage with point data types providing public access to x,y,z coordinates. `example2` shows how to use a different point type, which non-public coordinates. `example3` shows how to use the templated method inside an also templated descriptor.

We also provide a test case using the [Google Test Framework (GTest)](https://code.google.com/p/googletest/), which is automatically build if the package is eather found by Cmake or in the corresponding source directory.
You can invoke the testsuite with

```bash
./octree-test
```


## Attribution

If you use the implementation or ideas from the [corresponding paper](http://www.iai.uni-bonn.de/~behley/papers/behley2015icra.pdf) in your academic work, it would be nice if you cite the corresponding paper:

J. Behley, V. Steinhage, A.B. Cremers. *Efficient Radius Neighbor Search in Three-dimensional Point Clouds*, Proc. of the IEEE International Conference on Robotics and Automation (ICRA), 2015.

The BibTeX entry for the paper is::
    
    @conference{behley2015icra,
         author = {Jens Behley and Volker Steinhage and Armin B. Cremers},
          title = {{Efficient Radius Neighbor Seach in Three-dimensional Point Clouds}},
      booktitle = {Proc. of the IEEE International Conference on Robotics and Automation (ICRA)},
           year = {2015}
    }

## License


Copyright 2015 Jens Behley, University of Bonn.

This project is free software made available under the MIT License. For details see the LICENSE file.