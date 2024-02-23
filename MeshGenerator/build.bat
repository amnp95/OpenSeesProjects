@echo off
cd GKlib
REM Create directories
mkdir build
mkdir GKLIBRARY

REM Set environment variable
set GKLIBRARY=%CD%\GKLIBRARY

REM Change directory to build
cd build

REM Set GKLIBRARY environment variable to the full path of the GKLIBRARY directory

REM Assuming the second "cd build" is a mistake since we're already in "build" directory. 
REM If it's meant to go into another "build" directory inside the first, uncomment the next line:
REM cd build

REM Run cmake commands
cmake .. -DCMAKE_INSTALL_PREFIX=%GKLIBRARY%
cmake --build . --config Release --target GKlib
cmake --install . 
cd ..
mkdir  GKLIBRARY\include\win32
copy win32\adapt.h GKLIBRARY\include\win32
cd ..



set "GKLIBRARY=%CD%\GKlib\GKLIBRARY"
cd METIS
REM Set METISLIBRARY to the current directory followed by \METISLIBRARY
mkdir build
mkdir build\windows
mkdir build\xinclude
mkdir METISLIBRARY
set "METISLIBRARY=%CD%\METISLIBRARY"


REM Echo IDXWIDTH and REALWIDTH environment variables to metis.h
set "IDXWIDTH=#define IDXTYPEWIDTH 32"
set "REALWIDTH=#define REALTYPEWIDTH 32"
echo %IDXWIDTH% >> build\xinclude\metis.h
echo %REALWIDTH% >> build\xinclude\metis.h

REM Concatenate the contents of include\metis.h to build\xinclude\metis.h
type include\metis.h >> build\xinclude\metis.h

REM Copy include\CMakeLists.txt to build\xinclude
copy include\CMakeLists.txt build\xinclude

REM Change directory to build
cd build\windows

REM Run cmake commands with the previously set environment variables
cmake  -DCMAKE_INSTALL_PREFIX=%METISLIBRARY% -DGKLIB_PATH=%GKLIBRARY% ..\..
cmake --build . --config Release
cmake --install . --prefix %METISLIBRARY% 
cd ..\..\..

cd Partitioner
@REM cp Partition_Linux.c Partition.c
copy Partition_Linux.c Partition.c
mkdir build
cd build
cmake .. -DGKLIB_PATH=%GKLIBRARY% -DMETIS_PATH=%METISLIBRARY%
cmake --build . --config Release
cmake --install . 
cd ../..

python test.py




