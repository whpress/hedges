@echo off
REM Check if running in Visual Studio Developer Command Prompt
if "%VCINSTALLDIR%"=="" (
    echo Please run this batch file from the Visual Studio Developer Command Prompt.
    exit /b
)

REM Set the output directory for object files and executable to the 'build' directory
set OBJ_DIR=Build_Windows
set EXE_DIR=Build_Windows

REM Compile and link all source files, specifying the output directory for .obj files and the executable
cl /EHsc /Fo%OBJ_DIR%/ /Fe%EXE_DIR%/DNAcode_demo.exe DNAcode_demo.cpp DNAcode.cpp RSecc.cpp

REM Copy text file to EXE_DIR so test can be run from there
copy WizardOfOzInEsperanto.txt %EXE_DIR%\
