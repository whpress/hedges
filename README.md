#  [refactoring] HEDGES 


## TODO : numpy dev install in novel machine : bionic-2.7-numpy

test docker environnement, compilation and runtime with simple test
```
make test_docker
```

build and run test_programm
```
make build && make hedges_testprogramm
```


# HEDGES

A package for encoding and decoding arbitrary byte data to and from strands of DNA using a robust an error-correcting code (ECC).

### HEDGES Error-Correcting Code for DNA Storage Corrects Indels and Allows Sequence Constraints

**William H. Press, John A. Hawkins, Stephen Knox Jones Jr, Jeffrey M. Schaub, and Ilya J. Finkelstein**

*Proc Natl Acad Sci*. accepted for publication (June, 2020)

### Installation

The following instructions should work across platforms, except that installing virtualenv with apt-get is Ubuntu specific. For other platforms, install virtualenv appropriately if desired.

First, clone the repository to a local directory:

```
git clone https://github.com/whpress/hedges.git
```

Optionally, you can install into a virtual environment (recommended):

```
sudo apt-get install -y virtualenv
cd hedges
virtualenv envhedges
. envhedges/bin/activate
```

Now install required packages:

```
pip install numpy==1.13.3 && pip install -r requirements.txt && python setup.py install
```

### What is supplied
Supplied is not a single program, but a kit for variable user applications.  The kit consists of

1. C++ source code that compiles (in Linux or Windows) to the Python-includable module `NRpyDNAcode`.  Precompiled binaries are supplied for Python 2.7 in Linux and Windows, but recompilation may be necessary if these don't work.  This module implements the HEDGES "inner code" as described in the paper.

2.  C++ source code that compiles (in Linux or Windows) to the Python-includable module `NRpyRS`.  Precompiled binaries are supplied for Python 2.7 in Linux and Windows, but recompilation may be necessary if these don't work.  This module implements the Schifra Reed-Solomon Error Correcting Code Library.  See http://www.schifra.com  for details and license restrictions.  This module is not needed for the HEDGES inner code, but is needed only to implement the "outer code" as described in the paper.  Some users will instead want to utilize their own outer codes.
 
3.  Python program `print_module_test_files.py`, which verifies that the above modules can be loaded and prints their usage.  Most users will not need to use any of the routines in these files directly, but should instead use the Python functions in the following file:
 
4. Python program `test_program.py` .  This defines various user-level functions for implementing the HEDGES inner and Reed-Solomon outer codes as described in the paper.  The example inputs arbitrary bytes from the file `WizardOfOzInEsperanto.txt`, encodes a specified number of packets (each with 255 DNA strands), corrupts the strands with a specified level of random substitutions, insertions, and deletions, decodes the strands, and verifies the error correction.  To better validate the installation, the code rate and corruption level set by default are chosen to be stressful to HEDGES and is greater than that in an intended use case. 

### Testing and familiarization

Run the program `test_program.py` .  It should produce output comparable (but not identical) to the files `sample_linux_test_output.txt` and `sample_windows_test_output.txt`.  The output will not be identical, because different random numbers are used to create DNA errors in each run.

If the above works, then try varying some of the parameters.  In particular, you can change `coderatecode` to increase or decrease the code rate, the values `(srate,drate,irate)` to change the fraction of substitutions, deletions, and insertions generated for the test, and `totstrandlen`, the total strand length of the DNA (including left and right primers).  The many other parameters are either self-explanatory, or else described in the paper.  Most users will not initially need to change them.

### Recompiling the C++ modules

The modules are built using the Numerical Recipes C++ class library `nr3python.h` . This is included here and also freely available for unlimited distribution at http://numerical.recipes/nr3python.h .  Generally, you will not need to understand this library, but, if you are curious, a tutorial on its use is at http://numerical.recipes/nr3_python_tutorial.html .  You should also consult this tutorial if you have difficulty recompiling the modules.  Note that while other Numerical Recipes routines are copyright and require a license, no restricted routines are used in the two modules here supplied.

In Linux, go to the directory `LinuxC++Compile` containing the source code and run the script `compile_all.sh` .  Then copy the two files produced, `NRpyDNAcode.so` and `NRpyRS.so`, to the directory containing `test_program.py`.  The most common source of errors is the compiler's inability to find required Python and Numpy include and library files that are part of your Python installation.  Unfortunately, we can't help you with that.

In Windows, go to the directory `WindowsC++Compile` and fire up the Community Visual Studio 2019 solution `NRpyDNAcode.sln` .  This should build the two files (in the `x64\Release` directory) `NRpyDNAcode.pyd` and `NRpyRS.pyd` .  Copy these to the directory containing `test_program.py`.   If this doesn't work, and you need to build your the Windows modules from scratch, then keep these points in mind:  You want to compile to produce .dll files (not .exe files), and you want to then simply rename these to .pyd.  As in Linux, a common source of errors is the compiler's inability to find required Python and Numpy include and library files that are part of your Python installation.  You'll need to locate them and set appropriate include directories.

> Written with [StackEdit](https://stackedit.io/).
