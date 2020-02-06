# Temporary README for AccSeqV9

## What is this repo about

AccSeqV9 is a DNA sequence alignment tool which uses our novel seeding, exact string matching algorithm, and a new sequence alignment algorithm GACT (GACT is built according to papers from Darwin, which is using FPGA) algorithm to help with faster DNA sequencing. This repo contains a working version for single end DNA reads sequencing.

## Next Steps

- Accelerate the code using OpenACC to support GPU and other accelerators
- Run experiments and tune parameters to achieve better performance (processing time, accuracy and sensitivity metrics )
- Currently, a new repo is under development and it has significant changes to the voting algorithm in this version, which would result in better sensitivity and accuracy
- The GACT algorithm will be replaced by a new algorithm we developed which is still in development, and the GACT implementation in this code base is having bugs, so we are switching back to edlib right now.

## Usage
Compiling:
```
module load pgi
mkdir build
cd build
# Make sure CC=pgcc and CXX=pgc++ here and zlib-dev is correctly configured
cmake ../ -DCMAKE_BUILD_TYPE=Debug
# Going back to project root
cd ..
```
Running:
```
cd runtime
ln -s ../build/{accidx,accaln} ./
./accidx yeast.fa
./accaln yeast.fa yeast-1.fq > yeast-1.sam
```
The dataset in runtime folder is used as an example, feel free to apply your own dataset

## Known issues
This code could not be compiled with RelWithDebInfo because of a compiler bug (PGI compiler won't compile it in Release mode, but the same code compiles in GCC and PGI compiler in Debug mode)
This code is also having issues with cmake (weird bug reported but not replied), `find_package(ZLIB REQUIRED)` will fail in some systems (even if you have setup your system correctly), so please make sure zlib.h could be found in folders in CPATH, and libz.a could be fould in LIBRARY_PATH