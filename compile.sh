#!/bin/bash
# bash script will run as in a sandbox, so cd here won't change the pwd outside
# remember to `chmod +x compile.sh` for this file or use `sh compile.sh`
[[ -d build ]] || mkdir build
# rm build/CMakeCache.txt # if the directory is moved or reconstructed, reset the cmake cache before make
cd ./build
cmake ..

if [[ "$OSTYPE" == "msys" ]]; then
    # for Windows system
    cmake --build .
    # Note: we need to change the compile flags in /C-FLEX2D/FEM/CMakeLists.txt
    # executable will be under ./buld/FEM/Debug/main.ext
elif [[ "$OSTYPE" == "darwin"* ]]; then
    # for Mac OS
    make
    # executable will be under ./build/FEM/main
fi
