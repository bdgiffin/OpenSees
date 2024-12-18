#!/bin/bash

rm -rf build
mkdir build
cd build
/opt/homebrew/Cellar/conan@1/1.65.0/bin/conan install .. --build missing
cmake .. -DMUMPS_DIR=$PWD/../../mumps/build -DOPENMPI=TRUE -DSCALAPACK_LIBRARIES=/usr/local/Cellar/scalapack/2.2.0_1/lib/libscalapack.dylib
#cmake --build . --config Release --target OpenSees --parallel 4
#cmake --build . --config Release --target OpenSeesPy
