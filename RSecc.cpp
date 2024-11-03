#include "nr3b.h"
#include "ran.h"
#include "reed_solomon_schifra.h"

SchifraCode<255, 32> rs;
Ran rann;

VecUchar rsencode(VecUchar &message) {
	if (message.size() != 255) {
		throw("rsencode requires input array of size exactly 255");
	}
	VecUchar codetext = rs.encode(message);
	return codetext;
}

struct Rsdecode_out {
	VecUchar decoded;
	Int errs_detected;
	Int errs_corrected;
	Int err_number;
	bool recoverable;
};
Rsdecode_out rsdecode(VecUchar &received, VecInt &locations) {
	Int errs_detected, errs_corrected, err_code;
	bool recoverable;
	VecUchar decoded = rs.decode(received, locations, errs_detected, errs_corrected, err_code, recoverable);
	Rsdecode_out result = { decoded, errs_detected, errs_corrected, err_code, recoverable };
	return result;
}

struct Makeerasures_out {
	VecUchar codeword;
	VecInt location;
};
Makeerasures_out makeerasures(VecUchar &incodeword, Int nerase) {
	VecUchar codeword(incodeword); // so won't alter incodeword
	Int i, merase = 0, nn = codeword.size();
	Doub p;
	VecInt location(nerase);
	for (i = 0; i < nn; i++) {
		p = Doub(nerase - merase) / Doub(nn - i);
		if (rann.doub() < p) {
			location[merase++] = i;
			codeword[i] += 1; // wraparound; anything to change it is ok
		}
	}
	Makeerasures_out result = { codeword, location };
	return result;
}

VecUchar makeerrors(VecUchar &incodeword, Int nerror) {
	VecUchar codeword(incodeword); // so won't alter incodeword
	Int i, merror = 0, nn = codeword.size();
	Doub p;
	for (i = 0; i < nn; i++) {
		p = Doub(nerror-merror) / Doub(nn-i);
		if (rann.doub() < p) {
			codeword[i] += 1; // wraparound; anything to change it is ok
			++merror;
		}
	}
	return codeword;
}


