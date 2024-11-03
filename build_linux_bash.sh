#!/bin/bash

# Create the output directory if it doesn't exist
OBJ_DIR="Build_Linux"
EXE_DIR="Build_Linux"
mkdir -p $OBJ_DIR

# Compile and link all source files, specifying the output directory for .o files and the executable
g++ -std=c++11 -c -o $OBJ_DIR/DNAcode_demo.o DNAcode_demo.cpp
g++ -std=c++11 -c -o $OBJ_DIR/DNAcode.o DNAcode.cpp
g++ -std=c++11 -c -o $OBJ_DIR/RSecc.o RSecc.cpp

# Link the object files into the executable
g++ $OBJ_DIR/DNAcode_demo.o $OBJ_DIR/DNAcode.o $OBJ_DIR/RSecc.o -o $EXE_DIR/DNAcode_demo

# Copy the text file to the executable directory so tests can be run from there
cp WizardOfOzInEsperanto.txt $EXE_DIR/

echo "Build complete. Executable is in $EXE_DIR/DNAcode_demo."
