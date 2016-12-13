#!/bin/bash
module swap PrgEnv-intel PrgEnv-gnu
module swap gcc gcc/6.1.0
module load cmake3
module load boost
mkdir -p build/eos
cd build/eos
cmake ../..
make
