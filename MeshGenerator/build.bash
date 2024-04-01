cd GKlib
mkdir build
mkdir GKLIBRARY

export GKLIBRARY=$PWD/GKLIBRARY

cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$GKLIBRARY
cmake --build . --config Release --target GKlib
cmake --install .
cd ../..

cd METIS
mkdir build
mkdir build/xinclude
mkdir METISLIBRARY
mkdir METISLIBRARY/include
mkdir METISLIBRARY/lib

export METISLIBRARY=$PWD/METISLIBRARY


IDXWIDTH="#define IDXTYPEWIDTH 32"
REALWIDTH="#define REALTYPEWIDTH 32"
echo $IDXWIDTH >> build/xinclude/metis.h
echo $REALWIDTH >> build/xinclude/metis.h

cat include/metis.h >> build/xinclude/metis.h

cp include/CMakeLists.txt build/xinclude

cd build
cmake -DCMAKE_INSTALL_PREFIX=$METISLIBRARY -DGKLIB_PATH=$GKLIBRARY ..
cmake --build . --config Release --target metis

cp libmetis/libmetis.a "$METISLIBRARY/lib/libmetis.a"
cp xinclude/metis.h "$METISLIBRARY/include/metis.h"


cd ../..


cd Partitioner
cp Partition_Linux.c Partition.c
mkdir build
cd build
cmake .. -DGKLIB_PATH=$GKLIBRARY -DMETIS_PATH=$METISLIBRARY
cmake --build . --config Release
cmake --install .
cd ../..


