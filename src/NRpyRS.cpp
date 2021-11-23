#include "nr3python.h"
#include "ran.h"
#include "reed_solomon_schifra.h"

SchifraCode<255, 32> rs;
Ran ran;

static PyObject* rsencode(PyObject *self, PyObject *pyargs) {
	NRpyArgs args(pyargs);
	if (args.size() != 1) {
		NRpyException("rsencode takes 1 argument only");
		return NRpyObject(0); // formerly NULL
	}
	if (PyArray_TYPE(args[0]) != PyArray_UBYTE) {
		NRpyException("rsencode requires array with dtype=uint8 \n");
		return NRpyObject(0);
	}
	VecUchar message(args[0]);
	if (message.size() != 255) {
		NRpyException("rsencode requires input array of size exactly 255");
		return NRpyObject(0);
	}
	VecUchar codetext = rs.encode(message);
	return NRpyObject(codetext);
}

static PyObject* rsdecode(PyObject *self, PyObject *pyargs) {
	NRpyArgs args(pyargs);
	VecUchar received;
	VecInt locations;
	if (args.size() == 1) {
		received = VecUchar(args[0]);
		locations.resize(0);
	} else if (args.size() == 2) {
		received = VecUchar(args[0]);
		locations = VecInt(args[1]);
	} else {
		NRpyException("rsencode takes 1 or 2 arguments");
		return NRpyObject(0);
	}
	if (PyArray_TYPE(args[0]) != PyArray_UBYTE) {
		NRpyException("rsdecode requires array with dtype=uint8 \n");
		return NRpyObject(0);
	}
	Int errs_detected, errs_corrected, err_code;
	bool recoverable;
	VecUchar decoded = rs.decode(received, locations, errs_detected, errs_corrected, err_code, recoverable);
	return NRpyTuple(
		NRpyObject(decoded),
		NRpyObject(errs_detected),
		NRpyObject(errs_corrected),
		NRpyObject(err_code),
		NRpyObject(Int(recoverable)),
		NULL
	);
}

static PyObject* makeerasures(PyObject *self, PyObject *pyargs) {
	NRpyArgs args(pyargs);
	VecUchar incodeword(args[0]);
	VecUchar codeword(incodeword); // so won't alter incodeword
	Int nerase(NRpyInt(args[1]));
	Int i, merase = 0, nn = codeword.size();
	Doub p;
	VecInt location(nerase);
	for (i = 0; i < nn; i++) {
		p = Doub(nerase - merase) / Doub(nn - i);
		if (ran.doub() < p) {
			location[merase++] = i;
			codeword[i] += 1; // wraparound; anything to change it is ok
		}
	}
	return NRpyTuple(
		NRpyObject(codeword),
		NRpyObject(location),
		NULL
	);
}

static PyObject* makeerrors(PyObject *self, PyObject *pyargs) {
	NRpyArgs args(pyargs);
	VecUchar incodeword(args[0]);
	VecUchar codeword(incodeword); // so won't alter incodeword
	Int nerror(NRpyInt(args[1]));
	Int i, merror = 0, nn = codeword.size();
	Doub p;
	for (i = 0; i < nn; i++) {
		p = Doub(nerror-merror) / Doub(nn-i);
		if (ran.doub() < p) {
			codeword[i] += 1; // wraparound; anything to change it is ok
			++merror;
		}
	}
	return NRpyObject(codeword);
}


// standard boilerplate 
static PyMethodDef NRpyRS_methods[] = {
	{ "rsencode", rsencode, METH_VARARGS,
	"codetext = rsencode(uint8_array_length_255)" },
	{ "rsdecode", rsdecode, METH_VARARGS,
	"(decoded, errs_detected, errs_corrected, err_code, recoverable) = rsdecode(uint8_array_length_255[, VecInt erasure_locations])" },
	{ "makeerasures", makeerasures, METH_VARARGS,
	"(newcodetext,locations) = makeerasures(codetext,nerasures)" },
	{ "makeerrors", makeerrors, METH_VARARGS,
	"newcodetext = makeerrors(codetext,nerrors)" },
	{ NULL, NULL, 0, NULL }
};
PyMODINIT_FUNC initNRpyRS(void) {
	import_array();
	Py_InitModule("NRpyRS", NRpyRS_methods);
}

