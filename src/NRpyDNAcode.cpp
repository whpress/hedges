#include "nr3python.h"
#include "heapscheduler.h"
#include "ran.h"

//  this version 7 is version 6 with bug fixed in decode_c
//  this version 6 doesn't increment salt, but actually finds allowed output chars
//  this version 5 improves DNA constraints and does "fill" when codetext len is specified
//  this Version 4 adds DNA output constraints for GC balance and homopolymer runs
//  this Version 3 adds primers, check for coderate, and check for revcomp

Doub ThisVersion = 7.01;

#define GF4char Uchar // semantically ACGT
#define GF4word VecUchar  // semantically string of ACGT
#define GF4reg Ullong // semantically a compressed GF4word
#define Mbit Uchar // semantically VARIABLE NUMBER of plaintext bits
#define VecMbit VecUchar  // message bits unpacked to variable

// Globals: these used to be nicely contained in structs, but are here
// exposed for easy setting and inter-function communication

// user adjustable
Int NSALT = 24; // change salt after this many message bits (thus protecting them)
Int MAXSEQ = 2500; // maximum number of vbits in a message (one-time work in setcoderate() )
Int NSTAK = 110000;  // initial size of list of hypotheses
Int HLIMIT = 1000000; // limit on number of hypotheses tried before failure

// not normally user-adjustable
Int NPREV = 8; // number of hashed previous bits
Int NSEQBITS = 10; // number of hashed sequence number bits
Int HSALT = 24; // number of hashed bits of salt
Int LPRIMER = 0; // number of left-primer chars, set by findprimersalt()
Int RPRIMER = 0; // number of right-primer chars, set by findprimersalt()

Ullong prevmask((Ullong(1) << NPREV) - 1);
Ullong seqnomask((Ullong(1) << NSEQBITS) - 1);
Ullong saltmask((Ullong(1) << HSALT) - 1);
GF4word leftprimer, rightprimer;
VecUllong primersalt;
VecInt pattarr(MAXSEQ+2,1); // contains number of bits in each vbit: 0, 1, or 2
VecUchar pattrn(1,Uchar(1)); // initialize to rate 0.5 (pattnumber=3)
Int npattrn = 1, lastpattnumber=3;
Int VSALT = NSALT; // number of vbits corresponding to NSALT, updated  by setcoderate()  
Int NSP = VSALT + LPRIMER; // updated by setcoderate()

Int DNAWINDOW = 12; // window in which DNA constraints imposed
Int MAXGC = 8; // max GC in window
Int MINGC = 4; // min GC in window
Int MAXRUN = 4; // max length of homopolymers
GF4reg dnawinmask((Ullong(1) << 2 * DNAWINDOW) - 1);
GF4reg dnaoldmask((Ullong(1) << 2 * (DNAWINDOW-1)) - 1); // used to set oldest to "A"
GF4reg acgtacgt(0x1b1b1b1b1b1b1b1bllu); // "ACGTACGTACGTACGT" used for initialization

Uchar dnac_ok[4]; // sadly, another global
Int dnacallowed(GF4reg prev) {
	// returns the number of allowed ACGTs and puts them in dnac_ok
	if (DNAWINDOW <= 0) { dnac_ok[0] = 0; dnac_ok[1] = 1; dnac_ok[2] = 2; dnac_ok[3] = 3; return 4; }
	Int ans, gccount, last=prev & 3, nrun=1;
	bool isrun = false;
	Ullong reg;
	// get GCcount
	reg = prev & dnaoldmask;
	reg = (reg ^ (reg >> 1)) & 0x5555555555555555ull; // makes ones for GC, zeros for AT
	// popcount inline:
	reg -= ((reg >> 1) & 0x5555555555555555ull);
	reg = (reg & 0x3333333333333333ull) + (reg >> 2 & 0x3333333333333333ull);
	gccount = ((reg + (reg >> 4)) & 0xf0f0f0f0f0f0f0full) * 0x101010101010101ull >> 56; // the popcount
	// is there a run and, if so, of what
	reg = (prev >> 2);
	while ((reg & 3) == last) {
		++nrun;
		if (nrun >= MAXRUN) { isrun = true; break; }
		reg >>= 2;
	}
	// the horrible logic tree:
	if (gccount >= MAXGC) {
		ans = 2;
		dnac_ok[0] = 0; // A is ok
		dnac_ok[1] = 3; // T is ok
		if (isrun) {
			if (last == 0) {
				ans = 1;
				dnac_ok[0] = 3; // only T ok
			} else if (last == 3) {
				ans = 1;
				dnac_ok[0] = 0; // only A ok
			}
		}
	}
	else if (gccount <= MINGC) {
		ans = 2;
		dnac_ok[0] = 1; // C is ok
		dnac_ok[1] = 2; // G is ok
		if (isrun) {
			if (last == 1) {
				ans = 1;
				dnac_ok[0] = 2; // only G ok
			} else if (last == 2) {
				ans = 1;
				dnac_ok[0] = 1; // only C ok
			}
		}
	}
	else {  // no GC constraints
		ans = 4;
		dnac_ok[0] = 0; // A is ok
		dnac_ok[1] = 1; // C is ok
		dnac_ok[2] = 2; // G is ok
		dnac_ok[3] = 3; // T is ok
		if (isrun) {
			ans = 3;
			for (int i = last; i < 3; i++) dnac_ok[i] = dnac_ok[i + 1];
		}
	}
	return ans;
}

Ranhash ranhash;
inline Int digest(Ullong bits, Int seq, Ullong salt, Int mod) {
	return Int(ranhash.int64(
		((((Ullong(seq) & seqnomask) << NPREV) | bits) << HSALT) | salt
	) % mod);
}

static PyObject* getversion(PyObject *self, PyObject *pyargs) {
	NRpyArgs args(pyargs);
	return NRpyObject(ThisVersion);
}

static PyObject* getparams(PyObject *self, PyObject *pyargs) {
	NRpyArgs args(pyargs);
	return NRpyTuple(
		NRpyObject(NSALT),
		NRpyObject(MAXSEQ),
		NRpyObject(NSTAK),
		NRpyObject(HLIMIT),
		NULL);
}

static PyObject* restoreparams(PyObject *self, PyObject *pyargs) {
	NRpyArgs args(pyargs);
	NSALT = 24;
	MAXSEQ = 2500;
	NSTAK = 110000;
	HLIMIT = 1000000;
	return NRpyObject(Int(0));
}

static PyObject* setparams(PyObject *self, PyObject *pyargs) {
	NRpyArgs args(pyargs);
	if (args.size() != 4) {
		NRpyException("setparams takes exactly 4 arguments");
		return NRpyObject(Int(1));
	}
	NSALT = NRpyInt(args[0]);
	MAXSEQ = NRpyInt(args[1]);
	NSTAK = NRpyInt(args[2]);
	HLIMIT = NRpyInt(args[3]);
	return NRpyObject(Int(0));
}

static PyObject* getdnaconstraints(PyObject *self, PyObject *pyargs) {
	NRpyArgs args(pyargs);
	return NRpyTuple(
		NRpyObject(DNAWINDOW),
		NRpyObject(MAXGC),
		NRpyObject(MINGC),
		NRpyObject(MAXRUN),
		NULL);
}

static PyObject* restorednaconstraints(PyObject *self, PyObject *pyargs) {
	NRpyArgs args(pyargs);
	DNAWINDOW = 12;
	MAXGC = 8;
	MINGC = 4;
	MAXRUN = 4;
	dnawinmask = (Ullong(1) << 2 * DNAWINDOW) - 1;
	dnaoldmask = (Ullong(1) << 2 * (DNAWINDOW-1)) - 1;

	return NRpyObject(Int(0));
}

static PyObject* setdnaconstraints(PyObject *self, PyObject *pyargs) {
	NRpyArgs args(pyargs);
	if (args.size() != 4) {
		NRpyException("setdnaconstraints takes exactly 4 arguments");
		return NRpyObject(Int(1));
	}
	DNAWINDOW = NRpyInt(args[0]);
	MAXGC = NRpyInt(args[1]);
	MINGC = NRpyInt(args[2]);
	MAXRUN = NRpyInt(args[3]);
	dnawinmask = (Ullong(1) << 2 * DNAWINDOW) - 1;
	dnaoldmask = (Ullong(1) << 2 * (DNAWINDOW-1)) - 1;
	return NRpyObject(Int(0));
}

// these are the rewards and penalties applied at each position
Doub reward = -0.13;
Doub substitution = 1.;
Doub deletion = 1.;
Doub insertion = 1.;
Doub dither = 0.;

static PyObject* getscores(PyObject *self, PyObject *pyargs) {
	NRpyArgs args(pyargs);
	return NRpyTuple(
		NRpyObject(reward),
		NRpyObject(substitution),
		NRpyObject(deletion),
		NRpyObject(insertion),
		NRpyObject(dither),
		NULL);
}

static PyObject* restorescores(PyObject *self, PyObject *pyargs) {
	NRpyArgs args(pyargs);
	reward = -0.13;
	substitution = 1.;
	deletion = 1.;
	insertion = 1.;
	dither = 0.;
	return NRpyObject(Int(0));
}

static PyObject* setscores(PyObject *self, PyObject *pyargs) {
	NRpyArgs args(pyargs);
	if (args.size() != 5) {
		NRpyException("setscores takes exactly 5 arguments");
		return NRpyObject(Int(1));
	}
	reward = NRpyDoub(args[0]);
	substitution = NRpyDoub(args[1]);
	deletion = NRpyDoub(args[2]);
	insertion = NRpyDoub(args[3]);
	dither = NRpyDoub(args[4]);
	return NRpyObject(Int(0));
}

// more globals
Ran ran; // (11015);
GF4char* codetext_g; // set in decode, used by init_from_predecessor
Int codetextlen_g;   // ditto
Int nhypo = 0;
HeapScheduler<Doub,Int> heap;
Int errcode=0;
Int nfinal, nnstak;
Doub finalscore;
Int finaloffset, finalseq;

void findprimersalt(const char* leftpr, const char* rightpr) { // set salt to match a leftprimer
	Int regout, i, k, np = Int(strlen(leftpr)), mp = Int(strlen(rightpr));
	char ACGT[] = "ACGTacgt";
	VecInt ACGTvalue(256,0);
	LPRIMER = np;
	RPRIMER = mp;
	leftprimer.resize(np);
	primersalt.resize(np);
	rightprimer.resize(mp);
	for (i = 0; i < 8; i++) ACGTvalue[ACGT[i]] = i % 4;
	for (k = 0; k < np; k++) leftprimer[k] = ACGTvalue[leftpr[k]];
	for (k = 0; k < mp; k++) rightprimer[k] = ACGTvalue[rightpr[k]];
	for (k = 0; k < np; k++) {
		for (i = 0; i < 100; i++) { // try up to 100 times
			regout = digest(Ullong(0), k, Ullong(i), 4);
			if (regout == leftprimer[k]) {
				primersalt[k] = i;
				break;
			}
		}
	}
}

Int vbitlen(Int nmb) {  // how long is message in vbits?  (patarr must already be set)
	Int ksize, nn=0;
	for (ksize=0;;ksize++) {  // how many Mbits do we need?
		if (nn >= nmb) break;
		if (ksize >= MAXSEQ) NRpyException("vbitlen: MAXSEQ too small");
		nn += pattarr[ksize];
	}
	return ksize;
}

static PyObject* minstrandlen(PyObject *self, PyObject *pyargs) {
	NRpyArgs args(pyargs);
	Int nbytes = NRpyInt(args[0]);
	Int len = vbitlen(8 * nbytes) + RPRIMER;
	return NRpyObject(len);
}

// one more global below (hypostack)

static PyObject* hashint(PyObject *self, PyObject *pyargs) {
	NRpyArgs args(pyargs);
	if (args.size() != 1) {
		NRpyException("hashint takes exactly 1 argument");
		return NRpyObject(Int(1));
	}
	Int nn = NRpyInt(args[0]);
	Int hash = ranhash.int32(Ullong(nn));
	return NRpyObject(hash);
}

void setcoderate_C(Int pattnumber, const char* leftpr, const char* rightpr) { // some standard patterns
	findprimersalt(leftpr,rightpr);
	if (pattnumber == 1) { // rate 0.75
		pattrn.resize(2); pattrn[0] = 2; pattrn[1] = 1;
		reward = -0.035;
	}
	if (pattnumber == 2) { // rate 0.6
		pattrn.resize(5); pattrn[0] = 2;
		pattrn[1] = pattrn[2] = pattrn[3] = pattrn[4] = 1;
		reward = -0.082;
	}
	if (pattnumber == 3) { // rate 0.5
		pattrn.resize(1); pattrn[0] = 1;
		reward = -0.127;
	}
	if (pattnumber == 4) { // rate 0.333
		pattrn.resize(3); pattrn[0] = pattrn[1] = 1;  pattrn[2] = 0;
		reward = -0.229;
	}
	if (pattnumber == 5) { // rate 0.25
		pattrn.resize(2); pattrn[0] = 1; pattrn[1] = 0;
		reward = -0.265;
	}
	if (pattnumber == 6) { // rate 0.166
		pattrn.resize(3); pattrn[0] = 1; pattrn[1] = pattrn[2] = 0;
		reward = -0.324;
	}
	pattarr.assign(MAXSEQ + 2, 1);
	npattrn = pattrn.size();
	for (int i = 0; i < MAXSEQ; i++) pattarr[i] = (i < LPRIMER ? 0 : pattrn[i % npattrn]);
	VSALT = vbitlen(NSALT);
	NSP = VSALT + LPRIMER;
}

static PyObject* setcoderate(PyObject *self, PyObject *pyargs) {
	NRpyArgs args(pyargs);
	if (args.size() != 3) {
		NRpyException("setcoderate takes exactly 3 arguments");
		return NRpyObject(Int(1));
	}
	Int pattnumber = NRpyInt(args[0]);
	if (pattnumber < 1 || pattnumber > 6) {
		NRpyException("setcoderate arg must be in range 1 to 6");
		return NRpyObject(Int(1));
	}
	const char *leftpr = NRpyCharP(args[1]);
	const char *rightpr = NRpyCharP(args[2]);
	setcoderate_C(pattnumber, leftpr, rightpr);
	lastpattnumber = pattnumber;
	return NRpyObject(Int(0));
}

VecMbit unpackvbits(const char *message, Int n, Int len) {
	Int i, j, nmb=8*n, k, k1, ksize;
	Uchar bit;
	ksize = MAX(vbitlen(nmb), len-RPRIMER); // aim for codetext of length len if possible
	VecMbit ans(ksize,Uchar(0));
	i = j = 0;
	for (k = 0; k < ksize; k++) {
		for (k1 = 0; k1 < pattarr[k]; k1++) {
			bit = (i < n ? (message[i] >> (7-j++)) & 1 : 0);
			if (j == 8) { j = 0; ++i; }
			ans[k] = (ans[k] << 1) | bit;
		}
	}
	return ans;
}

VecUchar packvbits(VecMbit &vbits, Int nmessbits) {
	Int i,j,k,k1, ksize = vbits.size(), nn=0;
	Uchar bit;
	if (ksize > MAXSEQ) throw("packvbits: MAXSEQ too small");
	for (k = 0; k < ksize; k++) nn += pattarr[k]; // number of bits
	nn = MIN(nn, nmessbits); // no more than the specified number of bits
	nn = (nn + 7) / 8; // number of bytes
	VecUchar ans(nn,Uchar(0));
	i = j = 0;
	for (k = 0; k < ksize; k++) {
		for (k1 = pattarr[k] - 1; k1 >= 0; k1--) {
			bit = (vbits[k] >> k1) & 1;
			ans[i] = ans[i] | (bit << (7-j++));
			if (j == 8) {
				j = 0;
				if (++i == nn) break; 
			}
		}
		if (i == nn) break;
	}
	return ans;
}

Int bytepopcount (Uchar byte) {
  // not being used, but might someday!
  static const Uchar NIBBLE_LOOKUP [16] =
  {0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4};
  return NIBBLE_LOOKUP[byte & 0x0F] + NIBBLE_LOOKUP[byte >> 4];
}

GF4word encode_C(const char *message, Int n, Int len=0) { //dnac
	Int regout;
	GF4word vbits = unpackvbits(message, n, len);
	Int k=0, nbits, mod, nm = vbits.size(); // number of variable bits encoded
	if (nm > MAXSEQ) throw("encode: MAXSEQ too small");
	GF4word codetext(nm + RPRIMER);
	Mbit messagebit;
	Ullong prevbits = 0, salt = 0, newsalt = 0; 
	GF4reg prevcode = acgtacgt; // initialize with no runs and balanced cg
	for (k = 0; k < nm; k++) { // on decoding, k is called seq
		messagebit = vbits[k];
		nbits = pattarr[k];
		if (k < LPRIMER) {
			salt = primersalt[k];
		}
		else if (k < NSP) {
			salt = 0;
			newsalt = ((newsalt << 1) & saltmask) ^ messagebit;
		}
		else if (k == NSP) {
			salt = newsalt; // time to update the salt
		}
		mod = (k < LPRIMER ? 4 : dnacallowed(prevcode));
		regout = digest(prevbits, k, salt, mod);
		regout = (regout + Uchar(messagebit)) % mod;
		codetext[k] = (k < LPRIMER ? regout : dnac_ok[regout]);
		prevbits = ((prevbits << nbits) & prevmask) | messagebit; // variable number
		prevcode = ((prevcode << 2) | codetext[k]) & dnawinmask;
	}
	for (k = 0; k < RPRIMER; k++) {
		codetext[k + nm] = rightprimer[k];
	}
	return codetext;
}

GF4word encode_C(VecUchar &message, Int len=0) {
	return encode_C((char*)(&message[0]), Int(message.size()),len);
}
GF4word encode_C(char *message, Int len=0) {
	return encode_C(message, Int(strlen(message)),len);
}

static PyObject* encode(PyObject *self, PyObject *pyargs) {
	NRpyArgs args(pyargs);
	Int len = 0;
	if (args.size() > 2) {
		NRpyException("encode takes 1 or 2 arguments");
		return NRpyObject(0); // formerly NULL
	}
	if (PyArray_TYPE(args[0]) != PyArray_UBYTE) NRpyException("encode requires array with dtype=uint8 \n");
	VecUchar message(args[0]);
	if (args.size() > 1) len = NRpyInt(args[1]);
	VecUchar codetext = encode_C(message,len);
	return NRpyObject(codetext);
}

static PyObject* encodestring(PyObject *self, PyObject *pyargs) {
	NRpyArgs args(pyargs);
	GF4word empty(0);
	Int len = 0;
	if (args.size() > 2) {
		NRpyException("encodestring takes 1 or 2 arguments");
		return NRpyObject(0); // formerly NULL
	}
	const char *message = NRpyCharP(args[0]);
	if (args.size() > 1) len = NRpyInt(args[1]);
	VecUchar codetext = encode_C(message,len);
	return NRpyObject(codetext);
}

struct Hypothesis; // forward declaration for next line
NRvector<Hypothesis>* hypostackp; // pointed to hypostack by init_heap_and_stack()

struct Hypothesis {
	Int predi; // index of predecessor in hypostackp
	Int offset; // next char in message
	Int seq; // my position in the decoded message (0,1,...)
	Doub score; // my -logprob score before update
	Mbit messagebit; // last decoded up to now
	Ullong prevbits, salt, newsalt;
	GF4reg prevcode;
	
	Hypothesis() {}
	Hypothesis(int) {} // so that can cast from zero in NRvector constructor

	Int init_from_predecessor (Int pred, Mbit mbit, Int skew) {
		bool discrep;
		Int regout, mod;
		Doub mypenalty;
		Ullong mysalt;
		Hypothesis* hp = &(*hypostackp)[pred]; // temp pointer to predecessor
		predi = pred;
		messagebit = mbit; // variable number
		seq = hp->seq + 1;
		if (seq > MAXSEQ) throw("init_from_predecessor: MAXSEQ too small");
		Int nbits = pattarr[seq];
		prevbits = hp->prevbits;
		salt = hp->salt;
		if (seq < LPRIMER) {
			mysalt = primersalt[seq];
		}
		else if (seq < NSP) {
			mysalt = salt;
			newsalt = ((hp->newsalt << 1) & saltmask) ^ messagebit; // variable bits overlap, but that's ok with XOR
		}
		else if (seq == NSP) {
			mysalt = salt = hp->newsalt; // time to update the salt
		}
		else mysalt = salt;
		offset = hp->offset + 1 + skew;
		if (offset >= codetextlen_g) return 0; // i.e., false
		// calculate predicted message under this hypothesis
		prevcode = hp->prevcode;
		mod = (seq < LPRIMER ? 4 : dnacallowed(prevcode));
		regout = digest(prevbits, seq, mysalt, mod);
		regout = (regout + Uchar(messagebit)) % mod;
		regout = (seq < LPRIMER ? regout : dnac_ok[regout]);
		prevbits = ((hp->prevbits << nbits) & prevmask) | messagebit; // variable number
		prevcode = ((prevcode << 2) | regout) & dnawinmask;
		// compare to observed message and score
		if (skew < 0) { //deletion
			mypenalty = deletion;			
		} else {
			discrep = (regout == codetext_g[offset]);  // the only place where a check is possible!
			if (skew == 0) mypenalty = (discrep ? reward : substitution);
			else { // insertion
				mypenalty = insertion + (discrep ? reward : substitution);
			}
		}
		if (dither > 0.) mypenalty += dither * (2.*ran.doub() - 1.);
		score = hp->score + mypenalty;
		return 1; // i.e., true
	}
	void init_root() {
		predi = -1;
		offset = -1;
		seq = -1;
		messagebit = 0; // not really a message bit
		prevbits = 0;
		score = 0.;
		salt = 0;
		newsalt = 0;
		prevcode = acgtacgt;
	}	
	Int myid() { return Int(this - &((*hypostackp)[0])); } // my index in hypostack
};

// final global
NRvector<Hypothesis> hypostack;

void release() { // give back heap and hypostack memory
	heap.reinit();
	hypostack.resize(NSTAK, false);
}

void init_heap_and_stack() {
	hypostackp = &hypostack;
	if (nnstak < NSTAK) {
		nnstak = NSTAK;
		hypostack.resize(NSTAK, false);
	}
	hypostack[0].init_root();
	nhypo = 1;
	heap.rewind();
	heap.push(1.e10, 0);
}

void shoveltheheap(Int limit, Int nmessbits) {
	// given the heap, keep processing it until offset limit, hypothesis limit, or an error is reached
	Int qq, seq, nguess, qqmax=-1, ofmax= -1, seqmax = vbitlen(nmessbits);
	Uchar mbit;
	Doub currscore;
	Hypothesis *hp = NULL;
	errcode = 0;
	while (true) {
		currscore = heap.pop(qq);
		hp = &hypostack[qq];
		seq = hp->seq;
		if (seq > MAXSEQ) NRpyException("shoveltheheap: MAXSEQ too small");
		nguess = 1 << pattarr[seq + 1]; // i.e., 1, 2, or 4
		if (hp->offset > ofmax) { // keep track of farthest gotten to
			ofmax = hp->offset;
			qqmax = qq;
		}
		if (currscore > 1.e10) break; // heap is empty
		if (hp->offset >= limit-1) break; // errcode 0 (nominal success)
		if (nmessbits > 0 && seq >= seqmax-1) break; // ditto when no. of message bits specified
		if (nhypo > HLIMIT) {
			errcode = 2;
			nfinal = qqmax;
			return;
		}
		if (nhypo + 12 >= nnstak) {
			nnstak *= 2;
			hypostack.resize(nnstak, true);
			if (hypostack.size() != nnstak) NRpyException("resize of hypostack failed");
		}
		for (mbit = 0; mbit < nguess; mbit++) {
			if (hypostack[nhypo].init_from_predecessor(qq, mbit, 0)) { // substitution
				heap.push(hypostack[nhypo].score, nhypo);
				nhypo++;
			}
		}
		for (mbit = 0; mbit < nguess; mbit++) {
			if (hypostack[nhypo].init_from_predecessor(qq, mbit, -1)) { // deletion
				heap.push(hypostack[nhypo].score, nhypo);
				nhypo++;
			}
		}
		for (mbit = 0; mbit < nguess; mbit++) {
			if (hypostack[nhypo].init_from_predecessor(qq, mbit, 1)) { // insertion
				heap.push(hypostack[nhypo].score, nhypo);
				nhypo++;
			}
		}
	}
	nfinal = qq; // final position
}

VecMbit traceback() {
	Int k, kk=0, q=nfinal;
	while ((q = hypostack[q].predi) > 0) ++kk; // get length of chain
	VecMbit ans(kk + 1); // each with variable bits
	finalscore = hypostack[nfinal].score;
	finaloffset = hypostack[nfinal].offset;
	finalseq = hypostack[nfinal].seq;
	q = nfinal;
	k = kk;
	ans[k--] = hypostack[q].messagebit;
	while ((q = hypostack[q].predi) > 0) {
		ans[k] = hypostack[q].messagebit;
		--k;
	}
	return ans;
}

// global containers for fulldata
VecInt allseq;
VecInt allnhypo;
VecInt alloffset;
VecDoub allscore;
VecInt allpredi;
VecUchar allmessagebit;
VecInt allprevbits;
VecInt allsalt;
VecInt allnewsalt;

void traceback_fulldata(NRvector<Hypothesis> &hypostack) {
// TODO: questionable! messagebit might be 0, 1 or 2 bits.  how are you supposed to know?
// see packvbits()
	Int k, kk = 0, q = nfinal;
	while ((q = hypostack[q].predi) > 0) ++kk; // get length of chain
	finalscore = hypostack[nfinal].score;
	finaloffset = hypostack[nfinal].offset;
	finalseq = hypostack[nfinal].seq;
	allseq.resize(kk + 1);
	alloffset.resize(kk + 1);
	allscore.resize(kk + 1);
	allnhypo.resize(kk + 1);
	allpredi.resize(kk + 1);
	allmessagebit.resize(kk + 1);
	allprevbits.resize(kk + 1);
	allsalt.resize(kk + 1);
	allnewsalt.resize(kk + 1);
	finalscore = hypostack[nfinal].score;
	finaloffset = hypostack[nfinal].offset;
	finalseq = hypostack[nfinal].seq;
	q = nfinal;
	k = kk;
	allseq[k] = hypostack[q].seq;
	alloffset[k] = hypostack[q].offset;
	allscore[k] = hypostack[q].score;
	allnhypo[k] = q;
	allpredi[k] = hypostack[q].predi;
	allmessagebit[k] = hypostack[q].messagebit;
	allprevbits[k] = Int(hypostack[q].prevbits); // only returning 32 (or 31) bits of these
	allsalt[k] = Int(hypostack[q].salt);
	allnewsalt[k] = Int(hypostack[q].newsalt);
	--k;
	while ((q = hypostack[q].predi) > 0) {
		allseq[k] = hypostack[q].seq;
		alloffset[k] = hypostack[q].offset;
		allscore[k] = hypostack[q].score;
		allnhypo[k] = q;
		allpredi[k] = hypostack[q].predi;
		allmessagebit[k] = hypostack[q].messagebit;
		allprevbits[k] = Int(hypostack[q].prevbits); // only returning 32 (or 31) bits of these
		allsalt[k] = Int(hypostack[q].salt);
		allnewsalt[k] = Int(hypostack[q].newsalt);
		--k;
	}
}

static PyObject* releaseall(PyObject *self, PyObject *pyargs) {
	NRpyArgs args(pyargs);
	heap.reinit();
	hypostack.resize(NSTAK);
	allseq.resize(0);
	alloffset.resize(0);
	allscore.resize(0);
	allpredi.resize(0);
	allmessagebit.resize(0);
	allprevbits.resize(0);
	allsalt.resize(0);
	allnewsalt.resize(0);
	return NRpyObject(Int(0));
}

VecUchar decode_C(GF4word &codetext, Int nmessbits=0) {
	codetext_g = &codetext[0]; // set the pointer
	codetextlen_g = codetext.size();
	init_heap_and_stack();
	shoveltheheap(codetext.size(), nmessbits); // THIS WAS BUG: //last arg was nmessbits, but now always do whole codetext
	VecMbit trba = traceback();
	VecUchar pack = packvbits(trba, nmessbits); // truncate only at the end
	return pack;
}

void decode_fulldata_C(GF4word codetext) {
	codetext_g = &codetext[0]; // set the pointer
	codetextlen_g = codetext.size();
	init_heap_and_stack();
	shoveltheheap(codetext.size(), 0);
	traceback_fulldata(hypostack);
}

static PyObject* decode(PyObject *self, PyObject *pyargs) {
	NRpyArgs args(pyargs);
	Int nmessbits;
	if (args.size() == 1) {
		nmessbits = 0;
	}
	else if (args.size() == 2) {
		nmessbits = NRpyInt(args[1]);
	}
	else {
		NRpyException("decode takes 1 or 2 arguments only");
		return NRpyObject(0); // formerly NULL
	}
	if (PyArray_TYPE(args[0]) != PyArray_UBYTE) NRpyException("decode requires array with dtype=uint8 \n");
	GF4word codetext(args[0]);
	VecUchar plaintext = decode_C(codetext,nmessbits);
	return NRpyTuple(
		NRpyObject(errcode),
		NRpyObject(plaintext),
		NRpyObject(nhypo),
		NRpyObject(finalscore),
		NRpyObject(finaloffset),
		NRpyObject(finalseq),
		NULL);
}

static PyObject* decode_fulldata(PyObject *self, PyObject *pyargs) {
	NRpyArgs args(pyargs);
	if (args.size() != 1) {
		NRpyException("decode takes 1 argument only");
		return NRpyObject(0); // formerly NULL
	}
	if (PyArray_TYPE(args[0]) != PyArray_UBYTE) NRpyException("decode requires array with dtype=uint8 \n");
	GF4word codetext(args[0]);
	decode_fulldata_C(codetext);
	VecUchar t_allmessagebit(allmessagebit);
	VecInt t_allseq(allseq), t_alloffset(alloffset), t_allpredi(allpredi), t_allprevbits(allprevbits),
	  t_allsalt(allsalt), t_allnewsalt(allnewsalt), t_allnhypo(allnhypo);
	VecDoub t_allscore(allscore);
	return NRpyTuple(
		NRpyObject(errcode),
		NRpyObject(nhypo),
		NRpyObject(t_allmessagebit), // must return a temp, because Python gets control of its contents!
		NRpyObject(t_allseq),
		NRpyObject(t_alloffset),
		NRpyObject(t_allscore),
		NRpyObject(t_allnhypo),
		NRpyObject(t_allpredi),
		NRpyObject(t_allprevbits),
		NRpyObject(t_allsalt),
		NRpyObject(t_allnewsalt),
		NULL
	);
}

// in-place reverse complement for GF4word
void revcomp_C(GF4word &arr) {
	Int i, len = arr.size();
	Uchar TGCA[] = { 3,2,1,0 };
	for (i=0;i<len/2;i++) SWAP(arr[i],arr[len-1-i]);
	for (i=0;i<len;i++) arr[i] = (arr[i] > 3 ? arr[i] : TGCA[arr[i]]);
}

static PyObject* revcomp(PyObject *self, PyObject *pyargs) {
	NRpyArgs args(pyargs);
	if (args.size() != 1) {
		NRpyException("revcomp takes exactly 1 argument");
		return NRpyObject(0);
	}
	if (PyArray_TYPE(args[0]) != PyArray_UBYTE) NRpyException("revcomp requires array with dtype=uint8 \n");
	GF4word arr(args[0]);
	revcomp_C(arr);
	return NRpyObject(0);
}

VecInt gethowfar(Int hlimit, Int maxseq, GF4word &codetext, const char* leftpr, const char* rightpr) {
	VecInt ans(7,0); // pattern 0 is not defined
	Int HLIMIT_save = HLIMIT, MAXSEQ_save = MAXSEQ, pattno_save = lastpattnumber;
	HLIMIT = hlimit; // change the globals
	MAXSEQ = maxseq;
	VecUchar dc;
	for (int ipatt = 1; ipatt <= 6; ipatt++) {
		setcoderate_C(ipatt, leftpr, rightpr);
		dc = decode_C(codetext);
		ans[ipatt] = finaloffset;
	}
	HLIMIT = HLIMIT_save; // restore the globals
	MAXSEQ = MAXSEQ_save;
	lastpattnumber = pattno_save;
	setcoderate_C(lastpattnumber, leftpr, rightpr);
	return ans;
}

static PyObject* tryallcoderates(PyObject *self, PyObject *pyargs) {
	NRpyArgs args(pyargs);
	if (args.size() != 5) {
		NRpyException("tryallcoderates takes exactly 5 arguments");
		return NRpyObject(0);
	}
	Int hlimit = NRpyInt(args[0]);
	Int maxseq = NRpyInt(args[1]);
	if (PyArray_TYPE(args[2]) != PyArray_UBYTE) NRpyException("tryallcoderates requires array with dtype=uint8 \n");
	GF4word codetext(args[2]);
	const char *leftpr = NRpyCharP(args[3]);
	const char *rightpr = NRpyCharP(args[4]);
	VecInt maxoffsets = gethowfar(hlimit, maxseq, codetext, leftpr, rightpr);
	return NRpyObject(maxoffsets);
}

static PyObject* createerrors(PyObject *self, PyObject *pyargs) {
	NRpyArgs args(pyargs);
	if (args.size() != 4) {
		NRpyException("createerrors takes exactly 4 arguments");
		return NRpyObject(0); // formerly NULL
	}
	if (PyArray_TYPE(args[0]) != PyArray_UBYTE) NRpyException("createrrors requires array with dtype=uint8 \n");
	GF4word codetext(args[0]);
	Doub srate = NRpyDoub(args[1]); 
	Doub drate = NRpyDoub(args[2]); 
	Doub irate = NRpyDoub(args[3]); 
	Int n=0, nn=codetext.size(), k = 0;
	GF4word ans(2 * nn); // overkill
	while (n < nn) {
		if (ran.doub() < irate) { // insertion
			ans[k++] = ran.int32() % 4;
			continue;
		}
		if (ran.doub() < drate) { // deletion
			++n;
			continue;
		}
		if (ran.doub() < srate) { //substitution or errorfree
			ans[k++] = (codetext[n++] + (ran.int32() % 3) + 1) % 4;
		} else {
			ans[k++] = codetext[n++];
		}
	}
	ans.resize(k, true);
	return NRpyObject(ans);
}

Doub primerscore(const char *ain, GF4word &bin, Int binlen) {
	// returns penalty of match (large is bad)
	Int i, j;
	Doub dn, rt, dg;
	Doub mispen = 1., gappen = 1., skwpen = 1.;
	char ACGT[] = "ACGT";
	Int ia = Int(strlen(ain)), ib = binlen;
	MatDoub cost(ia + 1, ib + 1);
	cost[0][0] = 0.;
	for (i = 1; i <= ia; i++) cost[i][0] = cost[i - 1][0] + skwpen;
	for (i = 1; i <= ib; i++) cost[0][i] = cost[0][i - 1] + skwpen;
	for (i = 1; i <= ia; i++) for (j = 1; j <= ib; j++) {
		dn = cost[i - 1][j] + ((j == ib) ? skwpen : gappen);
		rt = cost[i][j - 1] + ((i == ia) ? skwpen : gappen);
		dg = cost[i - 1][j - 1] + ((ain[i - 1] == ACGT[bin[j - 1]]) ? -1. : mispen);
		cost[i][j] = MIN(MIN(dn, rt), dg);
	}
	return cost[ia][ib];
}

void makesense_C(const char* leftprimer, GF4word &codeword) {
	// reverse complement codeword (in place) if that makes leftprimer agree better
	Int len = Int(strlen(leftprimer));
	Doub lscore = primerscore(leftprimer, codeword, len);
	GF4word rcodeword(codeword);
	revcomp_C(rcodeword);
	Doub rscore = primerscore(leftprimer, rcodeword, len);
	if (rscore <= lscore) {
		codeword = rcodeword;
	}
}

static PyObject* makegoodsense(PyObject *self, PyObject *pyargs) {
	NRpyArgs args(pyargs);
	if (args.size() != 2) {
		NRpyException("makegoodsense takes exactly 2 arguments");
		return NRpyObject(0);
	}
	const char *leftprimer = NRpyCharP(args[0]);
	if (PyArray_TYPE(args[1]) != PyArray_UBYTE) NRpyException("makegoodsense requires array with dtype=uint8 \n");
	GF4word codeword(args[1]);
	GF4word newcodeword(codeword);
	makesense_C(leftprimer, newcodeword);
	return NRpyObject(newcodeword);
}


// standard boilerplate 
static PyMethodDef NRpyDNAcode_methods[] = {
	{ "getversion", getversion, METH_VARARGS,
	"version = getversion()\n get version number as a float" },
	{ "minstrandlen", minstrandlen, METH_VARARGS,
	"minstrandlength = minstrandlen(nbytes)\n get min length of DNA to encode nbytes (then use longer!)" },
	{ "getparams", getparams, METH_VARARGS,
	"(NSALT, MAXSEQ, NSTAK, HLIMIT) = getparams()\n get current values of parameters" },
	{ "restoreparams", restoreparams, METH_VARARGS,
	"restoreparams()\n restore parameters to default values\nNB! must be followed by setcoderate" },
	{ "setparams", setparams, METH_VARARGS,
	"errorcode = setparams(NSALT, MAXSEQ, NSTAK, HLIMIT)\n set new parameter values\nNB! must be followed by setcoderate" },
	{ "getdnaconstraints", getdnaconstraints, METH_VARARGS,
	"(DNAWINDOW, MAXGC, MINGC, MAXRUN) = getdnaconstraints()\n get current values of DNA constraints" },
	{ "restorednaconstraints", restorednaconstraints, METH_VARARGS,
	"restorednaconstraints()\n restore DNA constraints to default values" },
	{ "setdnaconstraints", setdnaconstraints, METH_VARARGS,
	"errorcode = setdnaconstraints(DNAWINDOW, MAXGC, MINGC, MAXRUN)\n set new DNA constraint values\nDNAWINDOW=0 for no constraints" },	
	{ "getscores", getscores, METH_VARARGS,
	"(reward,substitution,deletion,insertion,dither) = getscores()\n get current scoring parameters" },
	{ "restorescores", restorescores, METH_VARARGS,
	"restorescores()\n restore scoring parameters to default values" },
	{ "setscores", setscores, METH_VARARGS,
	"errorcode = setscores(reward,substitution,deletion,insertion,dither)\n set new scoring parameters" },
	{ "setcoderate", setcoderate, METH_VARARGS,
	"errorcode = setcoderate(number, leftprimer, rightprimer)\n\
	 set coderate to one of six values for number=1..6 (0.75, 0.6, 0.5, 0.333, 0.25, 0.166)" },
	{ "encode", encode, METH_VARARGS,
	"int8_dna_array = encode(int8_message_array [, strandlen])\n encode a message with runout to strandlen" },
	{ "encodestring", encodestring, METH_VARARGS,
	"int8_dna_array = encodestring(message_as_string)\n encode a message" },
	{ "decode", decode, METH_VARARGS,
	"(errcode, int8_message_array, nhypo, score, offset, seq) = decode(int8_dna_array[, nmessbits])\n\
	decode a message optionally limited to nmessbits message bits" },
	{ "tryallcoderates", tryallcoderates, METH_VARARGS,
	"maxoffsets = tryallcoderates(hlimit, maxseq, int8_dna_array, leftprimer, rightprimer)\n\
	maxoffsets[i] is maximum offset achieved in trying coderate i (in 1..6) limited by hlimit" },
	{ "decode_fulldata", decode_fulldata, METH_VARARGS,
	"(errcode,nhypo,messagebit,seq,offset,score,hypo,predi,prevbits,salt,newsalt) =\n\
	decode_fulldata(codetext[, nmessbits])" },
	{ "createerrors", createerrors, METH_VARARGS,
	"new_int8_dna_array = createerrors(int8_dna_array, subrate, delrate, insrate)\n\
	create Poisson random errors at specified rates"},
	{ "releaseall", releaseall, METH_VARARGS,
	"errcode = releaseall()\n release memory grabbed by decode_fulldata"},
	{ "revcomp", revcomp, METH_VARARGS,
	"revcomp(int8_dna_array)\n reverse-complement a dna array in place"},
	{ "makegoodsense", makegoodsense, METH_VARARGS,
	"new_int8_dna_array = makegoodsense(leftprimer, int8_dna_array)\n\
	return array or its reverse-complement, whichever agrees best with leftprimer"},
	{ "hashint", hashint, METH_VARARGS,
	"hashedint = hashint(int)\n hash an integer by same algorithm as used throughout"},
	{ NULL, NULL, 0, NULL }
};
PyMODINIT_FUNC initNRpyDNAcode(void) { // N.B. must rename to agree with module name
	import_array();
	Py_InitModule("NRpyDNAcode", NRpyDNAcode_methods);  // N.B. must rename first arg, not second
}

