#!/bin/bash

EXE_NAME=$1

COMPILER=$(root-config --cxx)
FLAGS=$(root-config --cflags --libs)

$COMPILER -g -O3 -Wall -Wextra -Wpedantic -o $EXE_NAME $EXE_NAME.cxx $FLAGS -lTreePlayer
