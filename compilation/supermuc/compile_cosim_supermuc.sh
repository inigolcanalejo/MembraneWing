module load cmake/3.10
module load gcc/7
module load python/2.7_intel

rm -rf build/
export CC=gcc
export CXX=g++
cmake -H"." -B"build" -DCO_SIM_IO_PYBIND="/hppfs/work/pn69ni/di73jef3/Softwares/pybind11" -DPYBIND11_PYTHON_VERSION=2.7
cmake --build "build" --target install

module unload cmake/3.10
module unload gcc/7
module unload python/2.7_intel