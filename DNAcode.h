// declarations for NR_DNAcode.cpp

// higher-level functions (used in test program)
// create a message and encode it in DNA
struct Createmesspacket_out { MatUchar packet; VecUchar plaintext; };
Createmesspacket_out createmesspacket(Int packno);
MatUchar protectmesspacket(const MatUchar& packetin);
MatUchar messtodna(const MatUchar& mpacket);

// introduce substitutions, deletions, and insertions as errors
MatUchar createerrors(const MatUchar& dnabag, Doub srate, Doub drate, Doub irate);

// decode the DNA back to error-corrected message
struct Dnatomess_out { MatUchar mpacket; MatUchar epacket; Int baddecodes; Int erasures; };
Dnatomess_out dnatomess(const MatUchar& dnapacket);
struct Correctmesspacket_out {
	MatUchar packet; Int tot_detect; Int tot_uncorrect; Int max_detect; Int max_uncorrect; Int toterrcodes;
};
Correctmesspacket_out correctmesspacket(const MatUchar& packetin, const MatUchar& epacket);
VecUchar extractplaintext(const MatUchar& cpacket);

// lower-level functions (can be used independently of the higher-level functions)
Doub getversion(); // returns HEDGES DNAcode version number
Int minstrandlen(Int nbytes); // returns minimum length of DNA strand implied by other parameters

struct Getparams_out { Int NSALT; Int MAXSEQ; Int NSTAK; Int HLIMIT; };
Getparams_out getparams(); // gets default (or current) settings of above parameters
void restoreparams(); // restore those parameters to defaults
void setparams(Int nsalt, Int maxseq, Int nstak, Int hlimit); //set those parameters

struct Getdnaconstraints_out { Int DNAWINDOW; Int MAXGC; Int MINGC; Int MAXRUN; };
Getdnaconstraints_out getdnaconstraints(); // gets default (or current) settings for DNA constraints
void restorednaconstraints(); // restore them to defaults
void setdnaconstraints(Int dnawindow, Int maxgc, Int mingc, Int maxrun); // set them

struct Getscores_out { Doub reward; Doub substitution; Doub deletion; Doub insertion; Doub dither; };
Getscores_out getscores(); // get default (or current) penalties (see PNAS paper)
void restorescores(); // restore them to default
void setscores(Doub reward_, Doub substitution_, Doub deletion_, Doub insertion_, Doub dither_); // set them

// IMPORTANT! Set coderate depending on the expected rate of errors:
void setcoderate(Int pattnumber, char* leftpr, char* rightpr); 

VecUchar encode(VecUchar message, Int len = 0); // HEDGES encode a message of bytes (as a vector)
VecUchar encodestring(const char* message, Int len = 0); // same, but message is a null-terminated string

void releaseall(); // release memory of the heap and elsewhere and start over

struct Decode_out {
	Int errcode; VecUchar plaintext; Int nhypo; Doub finalscore;
	Int finaloffset; Int finalseq;
};
Decode_out decode(VecUchar& codetext, Int nmessbits = 0); // decode DNA to plaintext, reporting stats

struct Decode_fulldata_out { // for debugging, return lots and lots of stuff
	Int errcode; Int nhypo; VecUchar t_allmessagebit; VecInt t_allseq; VecInt t_alloffset;
	VecDoub t_allscore; VecInt t_allnhypo; VecInt t_allpredi; VecInt t_allprevbits;
	VecInt t_allsalt; VecInt t_allnewsalt;
};
Decode_fulldata_out decode_fulldata(VecUchar& codetext);

void revcomp(VecUchar& arr); // utility: reverse complement a DNA vector
VecUchar createerrors(VecUchar codetext, Doub srate, Doub drate, Doub irate); // utility: create errors in DNA vector

// for debugging or diagnostic use, try decoding a DNA string with every
// available code rate, and return a vector of how far each one got before bombing
// from too many errors for that code rate
VecInt tryallcoderates(Int hlimit, Int maxseq, VecUchar codetext,
	const char* leftpr, const char* rightpr);

// utility: reverse complement codeword (in place) if that makes leftprimer agree better
// often the first step to do for sequenced DNA from a pool
VecUchar makegoodsense(const char* leftprimer, VecUchar& codeword);

