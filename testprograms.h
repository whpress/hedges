// main() test for NR_RS.cpp alone
// if this doesn't work for you, then the rest won't either
int main_test_RS() {
    Uint seed = static_cast<unsigned int>(time(0));
    Ran ran(seed);
    VecUchar message(255);
    for (int i = 0; i < 223; ++i) { // message length
        message[i] = ran.int8();
    }
    VecUchar encoded = rsencode(message);

    int nerror = 8;   // Number of errors to introduce
    int nerase = 16;   // Number of erasures to introduce
    VecUchar errored_message = makeerrors(encoded, nerror);
    Makeerasures_out erasures_result = makeerasures(errored_message, nerase);
    Rsdecode_out decoded_result = rsdecode(erasures_result.codeword, erasures_result.location);
    VecUchar decoded = decoded_result.decoded;
    // Display the results
    Int nbad = 0;
    for (int i = 0; i < message.size(); ++i) {
        printf("%3d: %2d %2d\n", i, message[i], decoded[i]);
        if (decoded[i] != message[i]) nbad++;
    }
    cout << "nbad: " << nbad << endl;
    cout << "Errors detected: " << decoded_result.errs_detected << endl;
    cout << "Errors corrected: " << decoded_result.errs_corrected << endl;
    cout << "Error code: " << decoded_result.err_number << endl;
    cout << "Message recoverable: " << (decoded_result.recoverable ? "Yes" : "No") << endl;
    return 0;
}

// from here on, the main test program

// DNAcode.cpp globals that can be set or accessed by main_test_all()

extern Int bytesperstrand;
extern Int strandsperpacketmessage;
extern Int messbytesperstrand;
extern Int messbytesperpacket;
extern Int strandIDbytes;
extern Int strandrunoutbytes;
extern Int totstrandlen;
extern Int leftlen;
extern Int rightlen;
extern Int strandlen;
// globals not normally user settable because assumed by Reed - Solomon outer code :
extern Int strandsperpacket;
extern Int strandsperpacketcheck;

int main_test_all() {
    Doub coderates_[] = { 0., 0.75, 0.6, 0.5, 1. / 3., 0.25, 1. / 6. }; // table of coderates 1..6
    VecDoub coderates(7, coderates_);

    // see the PNAS paper to understand the intended use of primers
    char leftprimer_s[] = "TCGAAGTCAGCGTGTATTGTATG"; // _s means "as a string"
    char rightprimer_s[] = "TAGTGAGTGCGATTAAGCGTGTT"; // for direct right appending (no reverse complement)

    // user-settable parameters for this test
    Int coderatecode = 3; // test this coderate in coderates table above
    Int npackets = 20; // number of packets (of 255 strands each) to generate and test
    Int hlimit = 1000000; // maximum size of decode heap, see paper

    // these lines are setting global variables in DNAcode.cpp
    totstrandlen = 300; // total length of DNA strand
    strandIDbytes = 2; // ID bytes each strand for packet and sequence number (is global)
    strandrunoutbytes = 2; // confirming bytes end of each strand (see paper)

    // sub,del,ins rates to simulate in this test (as multiple of our experimentally observed values):
    Doub ratefac = 1.5;
    Doub srate = ratefac * 0.0238;
    Doub drate = ratefac * 0.0082;
    Doub irate = ratefac * 0.0039;

    // set parameters for DNA constrants (normally not changed, except for no constraint)
    Int max_hpoly_run = 4; // max homopolymer length allowed (0 for no constraint)
    Int GC_window = 12; // window for GC count (0 for no constraint)
    Int max_GC = 8; // max GC allowed in window (0 for no constraint)
    Int min_GC = GC_window - max_GC;

    // get and reset parameters per above
    Getparams_out params = getparams();
    Int NSALT = params.NSALT;
    Int MAXSEQ = params.MAXSEQ;
    Int NSTAK = params.NSTAK;
    Int HLIMIT = params.HLIMIT;
    setparams(8 * strandIDbytes, MAXSEQ, NSTAK, hlimit); // change NSALT and HLIMIT
    setcoderate(coderatecode, leftprimer_s, rightprimer_s); // set code rate with left and right primers
    setdnaconstraints(GC_window, max_GC, min_GC, max_hpoly_run); // set DNA constraints (see paper)

    // set values of global variables derived from above
    leftlen = Int(strlen(leftprimer_s));
    rightlen = Int(strlen(rightprimer_s));
    strandlen = totstrandlen - leftlen - rightlen;
    strandsperpacketmessage = strandsperpacket - strandsperpacketcheck;
    bytesperstrand = Int(strandlen * coderates[coderatecode] / 4.);
    messbytesperstrand = bytesperstrand - strandIDbytes - strandrunoutbytes; // payload bytes per strand
    messbytesperpacket = strandsperpacket * messbytesperstrand; // payload bytes per packet of 255 strands
    
    printf("message bytes per packet = %d\n\n", messbytesperpacket);
 
    printf("for each packet, these statistics are shown in two groups:\n\
    1.1 HEDGES decode failures, 1.2 HEDGES bytes thus declared as erasures\n\
    1.3 R-S total errors detected in packet, 1.4 max errors detected in a single decode\n\
    2.1 R-S reported as initially-uncorrected-but-recoverable total, 2.2 same, but max in single decode\n\
    2.3 R-S total error codes; if zero, then R-S corrected all errors\n\
    2.4 Actual number of byte errors when compared to known plaintext input\n\n");
   
    Int badpackets = 0;
    VecInt Tots(8,Int(0));

    for (Int ipacket = 0; ipacket < npackets; ipacket++) {

        // encode
        Createmesspacket_out res_create = createmesspacket(ipacket); // plaintext to message packet
        MatUchar messpack(res_create.packet);
        VecUchar messplain(res_create.plaintext);
        MatUchar rspack = protectmesspacket(messpack); // Reed-Solomon protect the packet
        MatUchar dnapack = messtodna(rspack); // encode to strands of DNA containing payload messplain

        // simulate errors in DNA synthesis and sequencing
        MatUchar obspack = createerrors(dnapack, srate, drate, irate);

        // decode the strands
        Dnatomess_out res_dna = dnatomess(obspack);
        MatUchar dpacket(res_dna.mpacket), epacket(res_dna.epacket);
        Int baddecodes(res_dna.baddecodes), erasures(res_dna.erasures);

        Correctmesspacket_out res_cor = correctmesspacket(dpacket, epacket);
        MatUchar cpacket(res_cor.packet);
        Int tot_detect(res_cor.tot_detect), tot_uncorrect(res_cor.tot_uncorrect),
            max_detect(res_cor.max_detect), max_uncorrect(res_cor.max_uncorrect),
            toterrcodes(res_cor.toterrcodes);

        // check against ground truth
        VecUchar messcheck = extractplaintext(cpacket);
        Int badbytes = 0;
        for (Int i = 0; i < messplain.size(); i++)
            if (messplain[i] != messcheck[i]) ++badbytes;

        Int allstats_t[] = { baddecodes, erasures, tot_detect, max_detect,
            tot_uncorrect, max_uncorrect, toterrcodes, badbytes };
        Tots += VecInt(8, allstats_t);

        printf("%3d: (%3d %3d %3d %3d) (%3d %3d %3d %3d) ", ipacket, baddecodes, erasures,
            tot_detect, max_detect, tot_uncorrect, max_uncorrect, toterrcodes, badbytes);
        if (badbytes == 0) { printf("packet OK\n"); }
        else { printf("packet NOT ok\n"); badpackets += 1; }
    }
    if (badpackets == 0) { printf("ALL packets OK\n"); } else { printf("SOME packets NOT ok\n"); }
    printf("TOT: (%4d %4d %4d %4d) (%4d %4d %4d %4d)\n\n",
        Tots[0], Tots[1], Tots[2], Tots[3], Tots[4], Tots[5], Tots[6], Tots[7]);
    return 0;
}