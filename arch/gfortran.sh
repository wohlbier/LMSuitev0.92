#!/bin/bash
dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"/../

FC=gfortran

BUILD_TYPE=Debug
#BUILD_TYPE=RelWithDebInfo
#BUILD_TYPE=Release

FC=${FC} \
cmake \
  -DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
  $dir
