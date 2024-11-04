// search for "START HERE" (below) if you only want to peruse the higher-level functions!

// This DNAcode.cpp is a separate compilation unit with the internals of the HEDGES code.
// Once compiled, its functions (and some important globals) are available to user code.

//  this version 8 is the conversion to pure C++ (not Python interface)
//  this version 7 is version 6 with bug fixed in decode_c
//  this version 6 doesn't increment salt, but actually finds allowed output chars
//  this version 5 improves DNA constraints and does "fill" when codetext len is specified
//  this Version 4 adds DNA output constraints for GC balance and homopolymer runs
//  this Version 3 adds primers, check for coderate, and check for revcomp

#include "nr3b.h" // a version of Numerical Recipes class library
#include "ran.h" // Numerical Recipes code
#include "heapscheduler.h" // Numerical Recipes code
#include "RSecc.h" // Reed-Solomon declarations

// globals used by various functions; declared here, they must be set in main()
Int bytesperstrand;
Int strandsperpacketmessage;
Int messbytesperstrand;
Int messbytesperpacket;
Int strandIDbytes;
Int strandrunoutbytes;
Int totstrandlen;
Int leftlen;
Int rightlen;
Int strandlen;

Doub ThisVersion = 8.01;

#define GF4char Uchar // semantically ACGT has values 0,1,2,3
#define GF4word VecUchar  // semantically string of ACGT
#define GF4reg Ullong // semantically a compressed GF4word
#define Mbit Uchar // semantically VARIABLE NUMBER of plaintext bits
#define VecMbit VecUchar  // message bits unpacked to variable

// user adjustable (globals exposed by the linker)
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

// not adjustable because assumed by Reed-Solomon outer code :
Int strandsperpacket = 255;
Int strandsperpacketcheck = 32;

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

Uchar dnac_ok[4]; // global because used by multiple routines in this compilation unit
Int dnacallowed(GF4reg &prev) {
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

Ranhash ranhash; // in ran.h
inline Int digest(Ullong bits, Int seq, Ullong salt, Int mod) {
	return Int(ranhash.int64(
		((((Ullong(seq) & seqnomask) << NPREV) | bits) << HSALT) | salt
	) % mod);
}

Doub getversion() {
	return ThisVersion; // a global
}

struct Getparams_out {
	Int NSALT; Int MAXSEQ; Int NSTAK; Int HLIMIT;
};
Getparams_out getparams() {
	Getparams_out result = {NSALT, MAXSEQ, NSTAK, HLIMIT};
	return result;
}

void restoreparams() {
	NSALT = 24;
	MAXSEQ = 2500;
	NSTAK = 110000;
	HLIMIT = 1000000;
	return;
}

void setparams(Int nsalt, Int maxseq, Int nstak, Int hlimit) {
	NSALT = nsalt;
	MAXSEQ = maxseq;
	NSTAK = nstak;
	HLIMIT = hlimit;
	return;
}

struct Getdnaconstraints_out {
	Int DNAWINDOW; Int MAXGC; Int MINGC; Int MAXRUN;
};
Getdnaconstraints_out getdnaconstraints() {
	Getdnaconstraints_out result = { DNAWINDOW, MAXGC, MINGC, MAXRUN };
	return result;
}

void restorednaconstraints() {
	DNAWINDOW = 12;
	MAXGC = 8;
	MINGC = 4;
	MAXRUN = 4;
	dnawinmask = (Ullong(1) << 2 * DNAWINDOW) - 1;
	dnaoldmask = (Ullong(1) << 2 * (DNAWINDOW-1)) - 1;
	return;
}

void setdnaconstraints(Int dnawindow, Int maxgc, Int mingc, Int maxrun) {
	DNAWINDOW = dnawindow;
	MAXGC = maxgc;
	MINGC = mingc;
	MAXRUN = maxrun;
	dnawinmask = (Ullong(1) << 2 * DNAWINDOW) - 1;
	dnaoldmask = (Ullong(1) << 2 * (DNAWINDOW-1)) - 1;
	return;
}

// these are the rewards and penalties applied at each position
Doub reward = -0.13;
Doub substitution = 1.;
Doub deletion = 1.;
Doub insertion = 1.;
Doub dither = 0.;

struct Getscores_out {
	Doub reward; Doub substitution; Doub deletion; Doub insertion; Doub dither;
};
Getscores_out getscores() {
	Getscores_out result = {reward, substitution, deletion, insertion, dither};
	return result;
}

void restorescores() {
	reward = -0.13;
	substitution = 1.;
	deletion = 1.;
	insertion = 1.;
	dither = 0.;
	return;
}

void setscores(Doub reward_, Doub substitution_, Doub deletion_, Doub insertion_, Doub dither_) {
	reward = reward_;
	substitution = substitution_;
	deletion = deletion_;
	insertion = insertion_;
	dither = dither_;
	return;
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
		if (ksize >= MAXSEQ) THROW("vbitlen: MAXSEQ too small");
		nn += pattarr[ksize];
	}
	return ksize;
}

Int minstrandlen(Int nbytes) {
	Int len = vbitlen(8 * nbytes) + RPRIMER;
	return len;
}

// one more global below (hypostack)

Int hashint(Int nn) {
	Int hash = ranhash.int32(Ullong(nn));
	return hash;
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

void setcoderate(Int pattnumber, char* leftpr, char* rightpr) {
	if (pattnumber < 1 || pattnumber > 6) {
		THROW("setcoderate arg must be in range 1 to 6");
	}
	setcoderate_C(pattnumber, leftpr, rightpr);
	lastpattnumber = pattnumber;
	return;
} // WHP has some doubts about char* vs char vs const here

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
	if (ksize > MAXSEQ) THROW("packvbits: MAXSEQ too small");
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
	if (nm > MAXSEQ) THROW("encode: MAXSEQ too small");
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

GF4word encode_C(const VecUchar &message, Int len=0) {
	return encode_C((char*)(&message[0]), Int(message.size()), len);
}

GF4word encode_C(char *message, Int len=0) {
	return encode_C(message, Int(strlen(message)), len);
}

VecUchar encode(const VecUchar &message, Int len = 0) {
	VecUchar codetext = encode_C(message,len);
	return codetext;
}

VecUchar encodestring(const char *message, Int len = 0) {
	//GF4word empty(0);
	VecUchar codetext = encode_C(message,len);
	return codetext;
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
		if (seq > MAXSEQ) THROW("init_from_predecessor: MAXSEQ too small");
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
		if (seq > MAXSEQ) THROW("shoveltheheap: MAXSEQ too small");
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
			if (hypostack.size() != nnstak) THROW("resize of hypostack failed");
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

void releaseall() {
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
	return;
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

void decode_fulldata_C(GF4word &codetext) {
	codetext_g = &codetext[0]; // set the pointer
	codetextlen_g = codetext.size();
	init_heap_and_stack();
	shoveltheheap(codetext.size(), 0);
	traceback_fulldata(hypostack);
}

struct Decode_out {
	Int errcode; VecUchar plaintext; Int nhypo; Doub finalscore; Int finaloffset; Int finalseq;
};
Decode_out decode(GF4word &codetext, Int nmessbits = 0) {
	VecUchar plaintext = decode_C(codetext,nmessbits);
	Decode_out result = {errcode, plaintext, nhypo, finalscore, finaloffset, finalseq};
	return result;
}

struct Decode_fulldata_out {
	Int errcode; Int nhypo; VecUchar t_allmessagebit; VecInt t_allseq; VecInt t_alloffset;
	VecDoub t_allscore; VecInt t_allnhypo; VecInt t_allpredi; VecInt t_allprevbits;
	VecInt t_allsalt; VecInt t_allnewsalt;
};
Decode_fulldata_out decode_fulldata(GF4word &codetext) {
	decode_fulldata_C(codetext);
	VecUchar t_allmessagebit(allmessagebit);
	VecInt t_allseq(allseq), t_alloffset(alloffset), t_allpredi(allpredi), t_allprevbits(allprevbits),
	  t_allsalt(allsalt), t_allnewsalt(allnewsalt), t_allnhypo(allnhypo);
	VecDoub t_allscore(allscore);
	Decode_fulldata_out result = { errcode, nhypo, t_allmessagebit, t_allseq, t_alloffset,
		t_allscore, t_allnhypo, t_allpredi, t_allprevbits, t_allsalt, t_allnewsalt };
	return result;
}

// in-place reverse complement for GF4word
void revcomp_C(GF4word &arr) {
	Int i, len = arr.size();
	Uchar TGCA[] = { 3,2,1,0 };
	for (i=0;i<len/2;i++) SWAP(arr[i],arr[len-1-i]);
	for (i=0;i<len;i++) arr[i] = (arr[i] > 3 ? arr[i] : TGCA[arr[i]]);
}

void revcomp(GF4word &arr) {
	revcomp_C(arr);
	return;
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

VecInt tryallcoderates(Int hlimit, Int maxseq, GF4word codetext,
	const char *leftpr, const char *rightpr) {
	VecInt maxoffsets = gethowfar(hlimit, maxseq, codetext, leftpr, rightpr);
	return maxoffsets;
}

GF4word createerrors(GF4word &codetext, Doub srate, Doub drate, Doub irate) {
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
	return ans;
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

GF4word makegoodsense(const char* leftprimer, GF4word &codeword) {
	GF4word newcodeword(codeword);
	makesense_C(leftprimer, newcodeword);
	return newcodeword;
}

// START HERE if you only want to see the higher-level functions

class GetWizBytes { // make bytes in a file available for tests
public:
	VecUchar wizbytes;
	Int wizlen = 0;
	Int offset = 0;
	bool useFile;
	Ran ran;

	GetWizBytes(const char* filename) : offset(0), useFile(filename[0] != '\0') {
		if (useFile) {
			FILE* file = fopen(filename, "rb");
			if (file) {
				fseek(file, 0, SEEK_END);
				wizlen = ftell(file);
				fseek(file, 0, SEEK_SET);
				wizbytes.assign(wizlen, 0);
				fread(&wizbytes[0], 1, wizlen, file);
				fclose(file);
			}
			else {
				std::cerr << "Failed to open file: " << filename << std::endl;
				useFile = false; // Fallback to random if file open fails
			}
		}
	}
	VecUchar getBytes(Int n) {
		VecUchar bytes(n);
		if (offset + n > wizlen) { offset = 0; } // Reset offset if we reach the end
		if (useFile) {
			for (Int i = 0; i < n; ++i) { bytes[i] = wizbytes[offset + i]; }
		}
		else {
			for (Int i = 0; i < n; ++i) { bytes[i] = ran.int8(); }
		}
		offset += n;
		return bytes;
	}
};
GetWizBytes getwiz("WizardOfOzInEsperanto.txt"); // instantiate the stream

struct Createmesspacket_out { MatUchar packet; VecUchar plaintext; };
Createmesspacket_out createmesspacket(Int packno) {
	MatUchar packet(strandsperpacket, bytesperstrand, Uchar(0));
	VecUchar plaintext(strandsperpacketmessage * messbytesperstrand, Uchar(0));
	for (Int i = 0; i < strandsperpacket; ++i) {
		packet[i][0] = static_cast<Uchar>(packno);  // Assign packno to the first byte
		packet[i][1] = static_cast<Uchar>(i);       // Assign i to the second byte
		if (i < strandsperpacketmessage) {
			VecUchar ptext = getwiz.getBytes(messbytesperstrand); // assumes a global getwiz
			for (Int j = 0; j < messbytesperstrand; ++j) {
				packet[i][strandIDbytes + j] = ptext[j];
				plaintext[i * messbytesperstrand + j] = ptext[j];
			}
		}
	}
	return { packet, plaintext };
}

MatUchar protectmesspacket(const MatUchar& packetin) {
	MatUchar packet = packetin;
	VecUchar regin(strandsperpacket, Uchar(0));
	for (Int j = 0; j < messbytesperstrand; ++j) {
		for (Int i = 0; i < strandsperpacket; ++i) {
			regin[i] = packet[i][((j + i) % messbytesperstrand) + strandIDbytes];
		}
		VecUchar regout = rsencode(regin);
		for (Int i = 0; i < strandsperpacket; ++i) {
			packet[i][((j + i) % messbytesperstrand) + strandIDbytes] = regout[i];
		}
	}
	return packet;
}

MatUchar messtodna(const MatUchar& mpacket) {
	MatUchar dpacket(strandsperpacket, totstrandlen, Uchar(0));
	Uchar filler[] = { 0, 2, 1, 3, 0, 3, 2, 1, 2, 0, 3, 1, 3, 1, 2, 0, 2, 3, 1, 0, 3, 2, 1, 0, 1, 3 };
	for (Int i = 0; i < strandsperpacket; ++i) { // Loop through each strand in the packet
		VecUchar dna = encode(mpacket.row(i)); // how to get row of mpacket?
		if (dna.size() < totstrandlen) {
			VecUchar dnaleft(dna.slice(0, dna.size() - rightlen - 1));
			VecUchar dnaright(dna.slice(dna.size() - rightlen, dna.size() - 1));
			VecUchar fillerSegment(totstrandlen - dna.size());
			for (Int j = 0; j < fillerSegment.size(); ++j) { fillerSegment[j] = filler[j]; }
			dna = dnaleft.cat(fillerSegment).cat(dnaright);
		}
		dpacket.row(i) = dna;
	}
	return dpacket;
}

struct Dnatomess_out { MatUchar mpacket; MatUchar epacket; Int baddecodes; Int erasures; };
Dnatomess_out dnatomess(const MatUchar& dnapacket) {
	Int baddecodes = 0;
	Int erasures = 0;
	MatUchar mpacket(strandsperpacket, bytesperstrand, Uchar(0));
	MatUchar epacket(strandsperpacket, bytesperstrand, Uchar(1));
	for (Int i = 0; i < strandsperpacket; ++i) {
		VecUchar therow(dnapacket.row(i));
		Decode_out decoded = decode(therow, 8 * bytesperstrand);
		if (decoded.errcode > 0) {
			baddecodes += 1;
			erasures += MAX(0, messbytesperstrand - Int(decoded.plaintext.size()));
		}
		Int lenmin = MIN(Int(decoded.plaintext.size()), bytesperstrand);
		for (Int j = 0; j < lenmin; ++j) {
			mpacket[i][j] = decoded.plaintext[j];
			epacket[i][j] = 0; // Mark as non-erasure
		}
	}
	return { mpacket, epacket, baddecodes, erasures };
}

struct Correctmesspacket_out {
	MatUchar packet; Int tot_detect; Int tot_uncorrect; Int max_detect; Int max_uncorrect; Int toterrcodes;
};
Correctmesspacket_out correctmesspacket(const MatUchar& packetin, const MatUchar& epacket) {
	MatUchar packet = packetin;
	VecUchar regin(strandsperpacket, Uchar(0));
	VecUchar erase(strandsperpacket, Uchar(0));
	Int tot_detect = 0;
	Int tot_uncorrect = 0;
	Int max_detect = 0;
	Int max_uncorrect = 0;
	Int toterrcodes = 0;
	for (Int j = 0; j < messbytesperstrand; ++j) {
		Int num_erasures = 0;
		for (Int i = 0; i < strandsperpacket; ++i) {
			Int pos = ((j + i) % messbytesperstrand) + strandIDbytes;
			regin[i] = packet[i][pos];
			erase[i] = epacket[i][pos];
			num_erasures += erase[i];
		}
		VecInt locations(num_erasures);
		Int loc_index = 0;
		for (Int i = 0; i < strandsperpacket; ++i) {
			if (erase[i] != 0) {
				locations[loc_index++] = i;
			}
		}
		Rsdecode_out decoded = rsdecode(regin, locations);
		tot_detect += decoded.errs_detected;
		tot_uncorrect += MAX(0, decoded.errs_detected - decoded.errs_corrected);
		max_detect = MAX(max_detect, decoded.errs_detected);
		max_uncorrect = MAX(max_uncorrect, MAX(0, decoded.errs_detected - decoded.errs_corrected));
		toterrcodes += (decoded.err_number != 0 ? 1 : 0);
		// Update packet with decoded data
		for (Int i = 0; i < strandsperpacket; ++i) {
			Int pos = ((j + i) % messbytesperstrand) + strandIDbytes;
			packet[i][pos] = decoded.decoded[i];
		}
	}
	return { packet, tot_detect, tot_uncorrect, max_detect, max_uncorrect, toterrcodes };
}

VecUchar extractplaintext(const MatUchar& cpacket) {
	VecUchar plaintext(strandsperpacketmessage * messbytesperstrand, Uchar(0));
	for (Int i = 0; i < strandsperpacketmessage; ++i) {
		for (Int j = 0; j < messbytesperstrand; ++j) {
			plaintext[i * messbytesperstrand + j] = cpacket[i][strandIDbytes + j];
		}
	}
	return plaintext;
}

MatUchar createerrors(const MatUchar& dnabag, Doub srate, Doub drate, Doub irate) {
	Int nrows = dnabag.nrows();
	Int ncols = dnabag.ncols();
	MatUchar newbag(nrows, ncols, Uchar(0));
	for (Int i = 0; i < nrows; ++i) {
		VecUchar row_t(dnabag.row(i));
		VecUchar dna = createerrors(row_t, srate, drate, irate);
		Int lenmin = MIN(Int(dna.size()), ncols);
		for (Int j = 0; j < lenmin; ++j) {
			newbag[i][j] = dna[j];
		}
	}
	return newbag;
}



