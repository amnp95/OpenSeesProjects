@REM clean the builded files

@echo off
@REM  delete ./GKLib/build
rmdir /s /q GKLib\build
rmdir /s /q GKLib\GKLIBRARY

@REM  delete ./METIS/build
rmdir /s /q METIS\build
rmdir /s /q METIS\METISLIBRARY

@REM  delete ./Partitioner/build
rmdir /s /q Partitioner\build
