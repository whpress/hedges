/*
Test program for

HEDGES Error-Correcting Code for DNA Storage Corrects Indels and Allows Sequence Constraints
William H. Press, John A. Hawkins, Stephen Knox Jones Jr, Jeffrey M. Schaub, Ilya J. Finkelstein
Demonstration driver and installation validation program

We encode a specified number of packets from known plaintext, create DNA errors, then
decode the DNA and compare.
*/

#include "nr3b.h" // a version of the Numerical Recipes C++ class library
#include "ran.h" // Numerical Recipes random number routines

// declaration files for separate compilation units
#include "RSecc.h" // Reed-Solomon
#include "DNAcode.h" // HEDGES

// the test programs included here
#include "testprograms.h"

int main() {
	try {
		//main_test_RS(); // tests only the Reed-Solomon (sanity check)
		main_test_all(); // tests DNAcode (includes Reed-Solomon)
		return 0;
	}
	// simple error handling (see nr3b.h and change as desired) :
	catch (const exception& e) { 
		cerr << e.what() << std::endl;
		return 1;
	}
}

