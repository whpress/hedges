// declarations for NR_RS.cpp
struct Rsdecode_out { VecUchar decoded; Int errs_detected; Int errs_corrected; Int err_number; bool recoverable; };
struct Makeerasures_out { VecUchar codeword; VecInt location; };
VecUchar rsencode(VecUchar& message);
Rsdecode_out rsdecode(VecUchar& received, VecInt& locations);
Makeerasures_out makeerasures(VecUchar& incodeword, Int nerase);
VecUchar makeerrors(VecUchar& incodeword, Int nerror);

