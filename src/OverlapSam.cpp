/*By ANDREW FARRELL
 * OverlapSam.cpp
 * -------------------------------------------------- 
 * Assembles k-mers containing variation into contigs
 * that represent the variant sequence
 * --------------------------------------------------
 */

#include <fstream>
#include <iostream>
#include <math.h>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <sys/time.h>
#include <time.h>
#include <unordered_map>
#include <vector>


#include "Util.h"

using namespace std;
unordered_map<unsigned long, int> Mutations;
int HashSize = -1;
bool FullOut = false; // todo: get rid of this and use verbose per function
unordered_map<string, bool> DupCheck;

/***
 *
 * @param sequences
 * @param quals
 * @param Ap
 * @param A_qualp
 * @param Ai
 * @param overlap
 * @param index
 * @param minPercentPassed
 * @param PerfectMatch
 * @param MinOverlapPassed
 * @param Threads
 * @return
 */

/*
 * Ai = index in sequence vector
 * Compare sequence A to B, where B is the 9 subsequent sequences in the vector
 */
int Align3(vector<string> &sequences, vector<string> &quals, string Ap, string A_qualp, int Ai, int &overlap, int &index,
           float minPercentPassed, bool &PerfectMatch, int MinOverlapPassed, int Threads) {
    // int QualityOffset = 33; //=64; todo: UNUSED
    // int MinQual = 20; todo: UNUSED
    bool verbose = false;
    int bestScore = 0;
    //int NumReads = sequences.size(); todo: UNUSED

    // todo: these give overlapping windows? e.g. 1-10, 2-11, 3-12, etc.
    int start = Ai + 1;
    int end = start + 10;
    if (end > sequences.size()) {
        end = sequences.size();
    }

    if (FullOut) { cout << "staring alignment from " << start << " to " << end << endl; }

#pragma omp parallel for shared(Ap, A_qualp, index, overlap, bestScore) num_threads(Threads)
    for (int j = start; j < end; j++) {

        // Thread local variables
        int MinOverlap = MinOverlapPassed;
        float minPercent = minPercentPassed;
        int LocalBestScore = 0;
        int LocalIndex = -1;
        int LocalOverlap = 0;

        // todo: instead of this, pass into pragma as firstprivate
        string A = Ap;
        int A_length = A.length();
        string A_qual = A_qualp;

        string B;
        string B_qual;
        //#pragma omp critical
        {
            B = sequences[j]; // todo: reading from global here - indices are overlapping b/w threads so could be reading from same index at same time
            B_qual = quals[j];  // todo: reading from global here
        }
        int B_length = B.length();

        // Find smallest/largest window
        int window = -1;
        int longest = -1;
        bool A_smaller = true;
        if (B_length > A_length) {
            window = A_length;
            longest = B_length;
            A_smaller = false;
        } else {
            window = B_length;
            longest = A_length;
            A_smaller = true;
        }

        // todo: what is this
        int minMatchCount = window - (window * minPercent);
        int A_count = 0;
        int B_count = 0;

        // Difference b/w the two sequences
        int delta = longest - window;
        for (int i = 0; i <= delta; i++) {
            float score = 0;
            // First check bases within the window size
            for (int k = 0; k < window; k++) {
                // If we find a match that isn't N and has a quality score on both strands > 5, increment score
                if (A.c_str()[k + A_count] == B.c_str()[k + B_count]) {
                    if (B.c_str()[k + B_count] != 'N' && (int) A_qual.c_str()[k + A_count] > 5 &&
                        (int) B_qual.c_str()[k + B_count] > 5) {
                        score++;
                    }
                }
                // If we already have too few matches, abort this overlap
                if ((k - score) > minMatchCount) {
                    score = -1;
                    break;
                }
            }

            if (A_smaller) {
                A_count++;
            } else {
                B_count++;
            }

            if (verbose) {
                cout << "	Score = " << score << endl;
            }

            float percent = score / (window);
            if (percent >= minPercent) {
                if (FullOut) { cout << percent << " = " << score << " / " << window << endl; }

                if (LocalBestScore < score) {
                    LocalBestScore = score;
                    LocalIndex = j;
                    if (A_smaller) {
                        LocalOverlap = i * -1;
                    } else {
                        LocalOverlap = i;
                    }
                }
                if (score == window) {
                    PerfectMatch = true;
                    break;
                }
            }
        }

        // If we don't have a perfect match anywhere within the sliding window, keep trying to align outside of min length
        if (!PerfectMatch) {
            for (int i = window - 1; i >= MinOverlap; i--) {
                if (verbose) { cout << "i = " << i << endl; }

                float score = 0;
                for (int k = 0; k <= i; k++) {
                    if (verbose) {
                        cout << "	k = " << k << " so A = " << A_length - i + k << " /\\ B = " << 0 + k << endl;
                    }
                    if (verbose) {
                        cout << "		 A >> " << A.c_str()[A_length - i + k - 1] << "=" << B.c_str()[0 + k] << " << B"
                             << endl;
                    }

                    if (A.c_str()[A_length - i + k - 1] == B.c_str()[0 + k]) {
                        if (B.c_str()[0 + k] != 'N' && (int) A_qual.c_str()[A_length - i + k - 1] > 5 &&
                            (int) B_qual.c_str()[0 + k] > 5) {
                            score++;
                        }
                    }
                    if ((k - score) > minMatchCount) {
                        score = -1;
                        break;
                    }
                }
                if (verbose) { cout << "	Score = " << score << endl; }

                float percent = score / (k);
                if (percent >= minPercent) {
                    if (FullOut == true) {
                        cout << "percentsecond = " << percent << " = " << score << " / " << k << endl;
                    }
                    if (LocalBestScore < score) {
                        LocalBestScore = score;
                        LocalIndex = j;
                        LocalOverlap = i - A_length + 1;
                        if (score == i) {
                            break;
                        }
                    }
                }
            }

            for (int i = window - 1; i >= MinOverlap; i--) {
                if (verbose) { cout << "i = " << i << endl; }
                float score = 0;
                for (int k = 0; k <= i; k++) {
                    if (B.c_str()[B_length - i + k - 1] == A.c_str()[0 + k]) {
                        if (A.c_str()[0 + k] != 'N' && (int) A_qual.c_str()[0 + k] > 5 &&
                            (int) B_qual.c_str()[B_length - i + k - 1] > 5) {
                            score++;
                        }
                    }
                    if ((k - score) > minMatchCount) {
                        score = -1;
                        break;
                    }
                }

                if (verbose) { cout << "	Score = " << score << endl; }

                float percent = score / (k);
                if (percent >= minPercent) {
                    if (FullOut == true) {
                        cout << "percentthird = " << percent << " = " << score << " / " << k << endl;
                    }
                    if (LocalBestScore < score) {
                        LocalBestScore = score;
                        LocalIndex = j;
                        LocalOverlap = B_length - i - 1;
                        if (score == i) {
                            break;
                        }
                    }
                }
            }
        }
#pragma omp critical(updateCounts)
        {
            if (bestScore < LocalBestScore) {
                bestScore = LocalBestScore;
                index = LocalIndex;
                overlap = LocalOverlap;
            }
        }
    }
    return bestScore;
}

string ColapsContigs(string A, string B, int k, string A_qual, string &B_qual,
                     string &Ad, string &Bd, string &As, string &Bs) {
    bool verbose = false;
    if (verbose) {
        cout << "Combinding; \n" << A << endl << B << endl;
    }

    int Asize = A.length();
    int Bsize = B.length();

    int Aoffset = 0;
    int Boffset = 0;
    int window;
    string newString = "";
    string newQual = "";
    string newDepth = "";

    if (k > 0) {
        Aoffset = k;
    } else {
        Boffset = abs(k);
    }

    if (verbose) {
        cout << "K = " << k << " so Aofset = " << Aoffset
             << " and Boffset = " << Boffset << endl;
    }

    for (int i = 0; i < Asize + Bsize; i++) {
        char Abase = 'Z';
        char Bbase = 'Z';
        char A_qualual = '!';
        char B_qualual = '!';
        unsigned char Adep = 0;
        unsigned char Bdep = 0;

        if (((i - Aoffset) >= 0) && ((i - Aoffset) < A.length())) {
            Abase = A.c_str()[i - Aoffset];
            A_qualual = A_qual.c_str()[i - Aoffset];
            Adep = Ad.c_str()[i - Aoffset];
        } else {
            Abase = 'Z';
            A_qualual = '!';
            Adep = 0;
        }

        if (i - Boffset >= 0 && i - Boffset < B.length()) {
            Bbase = B.c_str()[i - Boffset];
            B_qualual = B_qual.c_str()[i - Boffset];
            Bdep = Bd.c_str()[i - Boffset];
        } else {
            Bbase = 'Z';
            B_qualual = '!';
            Bdep = 0;
        }

        if (verbose) {
            cout << "I = " << i << " Bi = " << i - Boffset << " Ai = " << i - Aoffset
                 << " thus " << Abase << "-" << Bbase << endl;
        }

        if (Abase == Bbase && Abase != 'Z') {
            newString += Abase;

            if (A_qualual >= B_qualual) {
                newQual += A_qualual;
            } else {
                newQual += B_qualual;
            }

            if ((int) Adep + (int) Bdep < 250) {
                newDepth += (Adep + Bdep);
            } else {
                newDepth += 250;
            }

        } else if (Abase == 'Z' && Bbase != 'Z') {
            newString += Bbase;
            newQual += B_qualual;
            newDepth += Bdep;
        } else if (Abase != 'Z' && Bbase == 'Z') {
            newString += Abase;
            newQual += A_qualual;
            newDepth += Adep;
        } else if (Abase != 'Z' && Bbase != 'Z') {

            if (Abase == 'N' && Bbase != 'N') {
                newString += Bbase;
                newQual += B_qualual;
                newDepth += Bdep;
            } else if (Abase != 'N' && Bbase == 'N') {
                newString += Abase;
                newQual += A_qualual;
                newDepth += Adep;
            } else if (A_qualual >= B_qualual) {
                newString += Abase;
                newQual += A_qualual;
                newDepth += Adep;
            } else {
                newString += Bbase;
                newQual += B_qualual;
                newDepth += Bdep;
            }

        } else if (Abase == 'Z' && Bbase == 'Z') {
            B_qual = newQual;
            Bd = newDepth;
            break;
        }
    }
    B_qual = newQual;
    Bd = newDepth;
    Bs += As;
    return newString;
}

// This needs to be re-written - assigning quality from passed pointer
// but actually returning string. Instead return a struct with both.
// todo
// Also make the if statements more clear and rename - not trimming ends, but trimming any bases that are N
string TrimNends(string S, string &qual) {
    bool base = false;
    string NewS = "";
    string NewQ = "";
    for (int i = S.size() - 1; i >= 0; i--) {
        if (base) {
            NewS = S.c_str()[i] + NewS;
            NewQ = qual.c_str()[i] + NewQ;
        } else if (S.c_str()[i] != 'A' && S.c_str()[i] != 'C' && S.c_str()[i] != 'G' && S.c_str()[i] != 'T') {}
        else {
            base = true;
            NewS = S.c_str()[i] + NewS;
            NewQ = qual.c_str()[i] + NewQ;
        }
    }
    qual = NewQ;
    return NewS;
}

string ReplaceLowQBase(string S, string qual, int min) {
    string NewS = "";
    for (int i = 0; i < S.size(); i++) {
        if (int(qual.c_str()[i]) - 33 < min) { NewS = NewS + 'N'; }
        else { NewS = NewS + S.c_str()[i]; }
    }
    return NewS;
}

string TrimKends(string S, string &qual, int TrimLen) {
    string NewS = "";
    string NewQ = "";
    if (S.size() - TrimLen - TrimLen > 0) {
        for (int i = TrimLen; i < S.size() - TrimLen; i++) {
            //		cout << "i = " << i << endl;
            NewS = NewS + S.c_str()[i];
            NewQ = NewQ + qual.c_str()[i];
        }
        //	cout << "TRIM CEHCK\n" << S << "\n" << "               " << NewS << "\n" << qual << "\n               " << NewQ << endl;
    } else {
        cout << "Warning read is shorter than the clipping length, just returning the read" << endl;
    }
    qual = NewQ;
    return NewS;
}

string TrimLowCoverageEnds(string S, string &quals, string &depth, int cutoff) {
    bool base = false;
    string NewS = "";
    string NewD = "";
    string NewQ = "";

    for (int i = S.size() - 1; i >= 0; i--) {
        if (base) {
            NewS = S.c_str()[i] + NewS;
            NewD = depth.c_str()[i] + NewD;
            NewQ = quals.c_str()[i] + NewQ;
        } else if ((int) depth.c_str()[i] > cutoff) {
            base = true;
            NewS = S.c_str()[i] + NewS;
            NewD = depth.c_str()[i] + NewD;
            NewQ = quals.c_str()[i] + NewQ;
        }
    }

    if (NewS.size() > 1) {
        S = NewS;
        depth = NewD;
        quals = NewQ;
        base = false;
        NewS = "";
        NewD = "";
        NewQ = "";

        for (int i = 0; i < S.size(); i++) {
            if (base) {
                NewS = NewS + S.c_str()[i];
                NewD = NewD + depth.c_str()[i];
                NewQ = NewQ + quals.c_str()[i];
            } else if ((int) depth.c_str()[i] > cutoff) {
                base = true;
                NewS = NewS + S.c_str()[i];
                NewD = NewD + depth.c_str()[i];
                NewQ = NewQ + quals.c_str()[i];
            }
        }
    }
    depth = NewD;
    quals = NewQ;
    return NewS;
}

string AdjustBases(string sequence, string qual) {
    int MinQ = 10;
    int QualOffset = 32;
    string NewString = "";
    for (int i = 0; i < sequence.length(); i++) {
        if (qual.c_str()[i] - QualOffset < MinQ) {
            NewString += 'N';
        } else {
            NewString += sequence.c_str()[i];
        }
    }

    if (NewString != sequence) {
        return NewString;
    }
}

bool replace(std::string &str, const std::string &from, const std::string &to) {
    size_t start_pos = str.find(from);
    if (start_pos == std::string::npos) {
        return false;
    }
    str.replace(start_pos, from.length(), to);
    return true;
}

bool validateFASTQD(string &L1, string &L2, string &L3, string &L4, string &L5,
                    string &L6) {
    if (L1.c_str()[0] != '@') {
        cout << "error header problems - " << L1 << endl;
        return false;
    }
    if (L2.size() != L4.size()) {
        cout << "error sequence and qual	problems - \n" << L2 << endl
             << L4 << endl;
        return false;
    }
    vector<string> temp = Util::Split(L6, ' ');
    if (temp.size() != L2.size()) {
        cout << "error counts	problems - " << L2.size() << " != " << temp.size()
             << endl;
        return false;
    }
    return true;
}

bool IsBitSet(int num, int bit) { return 1 == ((num >> bit) & 1); }

int GetReadOrientation(int flag) {
    bool is_set = IsBitSet(flag, 4);
    cout << "flag is " << flag << endl;
    cout << "orientation is " << is_set << endl;
    return is_set;
}

string FlipStrands(string strand) {
    string NewStrand = "";
    for (int i = 0; i < strand.size(); i++) {
        if (strand.c_str()[i] == '+')
            NewStrand += "-";
        else if (strand.c_str()[i] == '-')
            NewStrand += "+";
        else if (strand.c_str()[i] == '.')
            NewStrand += ".";
    }
    return NewStrand;
}

void compresStrand(string S, int &F, int &R) {
    for (int i = 0; i < S.size(); i++) {
        if (S.c_str()[i] == '+')
            F++;
        else if (S.c_str()[i] == '-')
            R++;
    }
    return;
}

// todo: get rid of global variables and pass them as arguments
// This method splits each hash by Ns and then checks if the hash is in the hash list
int CountHashes(string seq) {
    int count = 0;
    for (int i = 0; i < seq.size() - HashSize; i++) {
        string hash = seq.substr(i, HashSize);
        size_t found = hash.find("N");
        if (found == string::npos) {
            if (Mutations.count(Util::HashToLong(hash)) > 0) { count++; } // todo: is count ever >1?
        }
    }
    return count;
}

int NumLowQbases(string qual, int min) {
    int count = 0;
    for (int i = 0; i < qual.size(); i++) {
        if (int(qual.c_str()[i]) - 33 < min)
            count++;
    }
    return count;
}

/***
 *
 * @param argc 10
 * @param argv 1) BAM, 2) MinPercent, 3) MinOverlap, 4) MinCoverage, 5) SamFile Output??, 6) NodeStub, 7) LCcutoff, 8) HashPath, 9) Threads
 * @return 0
 */
int main(int argc, char *argv[]) {
    // Parse arguments
    float MinPercent = stof(argv[2]);
    long int MinOverlap = strtol(argv[3], nullptr, 0);
    long int MinCoverage = strtol(argv[4], nullptr, 0);
    long int LCcutoff = strtol(argv[7], nullptr, 0);
    string HashPath = argv[8];
    long int Threads = strtol(argv[9], nullptr, 0);

    // Check correct number or arguments
    cout << "you gave " << argc << " Arguments" << endl;
    if (argc != 10) {
        cout << "ERROR, incorrect number of arguments \n "
                "Call is: SAM, MinPercent, MinOverlap, MinCoverage, ReportStub, NodeStub, LCcutoff, Threads"
             << endl
             << "You Gave \n"
             << "File = " << argv[1]
             << "\n	MinPercent = " << argv[2]
             << "\n	MinOverlap = " << argv[3]
             << "\n	MinCoverage = " << argv[4]
             << "\n	FileStub = " << argv[5]
             << "\n	NodeStub = " << argv[6]
             << "\n	MinCov = " << argv[7]
             << "\n	HashPath = " << argv[8]
             << "\n	Threads = " << argv[9]
             << endl;
        return 0;
    }
    cout << "Correct number of arguments \n "
            "Call is: SAM, MinPercent, MinOverlap, MinCoverage, ReportStub, NodeStub LCcutoff Threads"
         << endl
         << "You Gave \n"
         << "File = " << argv[1]
         << "\n	MinPercent = " << argv[2]
         << "\n MinOverlap = " << argv[3]
         << "\n MinCoverage = " << argv[4]
         << "\n FileStub = " << argv[5]
         << "\n NodeStub = " << argv[6]
         << "\n MinCovTrimCov = " << argv[7]
         << "\n HashPath = " << argv[8]
         << "\n Threads = " << argv[9]
         << endl;

    // Open the file
    ifstream bam_file;
    bam_file.open(argv[1]);
    if (bam_file.is_open()) {
        cout << "File open - " << argv[1] << endl;
    } else {
        cout << "Error, Proband could not be opened";
        return 0;
    }


    ofstream report;
    string FirstPassFile = argv[1];
    std::stringstream ss;
    ss << argv[5] << ".fastq";
    FirstPassFile = ss.str();   // This is name of output file (will eventually be .sam.fastqd)
    report.open(FirstPassFile.c_str());
    cout << "min percent check = " << MinPercent << endl;
    if (report.is_open()) {
    } else {
        cout << "ERROR, Mut-Output file could not be opened - " << FirstPassFile << endl;
        return 0;
    }

    // todo: is this fastqd that gets output?
    ofstream Depreport;
    FirstPassFile += "d";
    Depreport.open(FirstPassFile.c_str());  // This is output stream for .sam.fastqd

    if (report.is_open()) {
    } else {
        cout << "ERROR, Mut-Output file could not be opened - " << FirstPassFile << endl;
        return 0;
    }

    string line;
    std::vector<string> sequences;
    std::vector<string> qual;
    std::vector<string> depth;
    std::vector<string> strand;
    std::vector<string> Unsequences;
    std::vector<string> Unqual;
    std::vector<string> Undepth;
    std::vector<string> Unstrand;

    // Correspond with fastqd format?
    int lines = -1;
    string L1;
    string L2;
    string L3;
    string L4;
    string L5;
    string L6;

    unsigned long LongHash;
    int Rejects = 0;
    string Fastqd = argv[1]; // todo: unused and misleading
    cout << "Reading in SAM \n";
    int counter = 0;
    int unalignedCounter = 0;
    int goodreads = 0;
    int lowMapQual = 0;
    int other = 0;

    cout << "reading in hash list" << endl;
    ifstream MutHashFile;
    MutHashFile.open(HashPath);
    if (MutHashFile.is_open()) {
        cout << "Parent File open - " << argv[1] << endl;
    } else {
        cout << "Error, ParentHashFile could not be opened";
        return 0;
    }


    // e.g. MutHashFile= $ProbandGenerator".k"$K"_c"$MutantMinCov".HashList
    // Populate Mutations lookup from Jellyfish hash list
    while (getline(MutHashFile, L1)) {

        // Split line by spaces
        vector<string> temp = Util::Split(L1, ' ');

        unsigned long b = Util::HashToLong(temp[0]);
        unsigned long revb = Util::HashToLong(Util::RevComp(temp[0]));
        Mutations.insert(pair<unsigned long, int>(b, 0));
        Mutations.insert(pair<unsigned long, int>(revb, 0));
        HashSize = temp[0].size(); // todo: this is only populated/used when last line assigned - is every hash size the same?
    }

    if (HashSize == -1) {
        cout << "ERROR Hash Size could not be determined by the HashFile" << endl;
        return -1;
    }
    cout << "HashSize = " << HashSize << endl;

    // Iterate through bam file
    while (getline(bam_file, L1)) {
        counter++;
        if (counter % 100 == 1) {
            cout << "Read in " << counter << " reads, with " << goodreads << " aligned reads, " << unalignedCounter
                 << " unaligned reads, " << lowMapQual << "low map qual and " << other << " other with rejected "
                 << Rejects
                 << " reads\r";
        }
        vector<string> temp = Util::Split(L1, '\t');    // Temporary string vector
        temp[9] = ReplaceLowQBase(temp[9], temp[10], 10);
        //temp[9] = TrimKends(temp[9], temp[10], 15);
        int ReadSize = temp[10].size();
        bool b[16];
        int v = atoi(temp[1].c_str()); // bitwise FLAG field

        // Generate a binary representation of FLAG
        // See samtools format for each bit meaning
        for (int j = 0; j < 16; ++j) {
            b[j] = 0 != (v & (1 << j));
        }
        //if (DupCheck.count(temp[9]) > 0)
        //{
        //	cout << "skipping exact match sequence " << L1 << endl;
        //}
        //else
        // todo: get rid of these brackets - confusing
        //{
            int lowq = NumLowQbases(temp[10], 20);
            int length = temp[10].length();
            //DupCheck[temp[9]] == true;

            // Reject if secondary alignment (b[8]) or supplementary alignment (b[11]) or PCR or optical duplicate (b[10])
            // or if the length of the seq is < 50 or if 33% of the bases are low quality
            if (b[8] or b[11] or b[10] or temp[9].length() < 50 or ((double) lowq / (double) length > 0.33)) {
                Rejects++;
            }
            // Unmapped bit is b[2]
            if (b[2])
                unalignedCounter++;
            else if (atoi(temp[4].c_str()) <
                     5) // temp[4] is MAPQ todo: make const var at top for mapq cutoff instead of hardcoding here
                lowMapQual++;
            else
                other++;
            string L4 = temp[10]; // QUAL field
            string L2 = temp[9]; // SEQ field

            // Remove non-ACGT bases from sequence (and their corresponding QUAL digit)
            // NOTE: this updates L2 AND L4 (even though doesn't return L4)
            L2 = TrimNends(L2, L4);

            int hashes = CountHashes(L2);
            //cout << "poorly mapped read " << L1 << endl;
            //cout << "with Hash = " << hashes << endl;
            // Only keep a read if the trimmed sequence is > 60% of the original sequence
            // todo: make 0.6 a const var at top
            if ((double) L2.size() / (double) ReadSize > .6) {
                ReadSize = L2.size();
                lines++;
                Unsequences.push_back(L2);
                Unqual.push_back(L4);
                if (hashes > 0) {
                    if (b[0] == 0) {
                        Unstrand.push_back(".");
                    } else if (b[4] == 0) {
                        Unstrand.push_back("+");
                    } else if (b[4] == 1) {
                        Unstrand.push_back("-");
                    }
                } else
                    Unstrand.push_back(".");

                string depths = "";
                unsigned char C = 1; // todo: this is misleading, don't really need a variable

                for (int i = 0; i < L2.length(); i++) {
                    depths += C;
                }

                Undepth.push_back(depths);

            } else {
                Rejects++;
            }

            // todo: I don't think this is ever executed since the if statement above is always true
//        } else {
//            goodreads++;
//            //cout << "good alignment" << endl;
//            string L4 = temp[10];
//            string L2 = temp[9];
//            L2 = TrimNends(L2, L4);
//            int hashes = CountHashes(L2);
//
//            //cout << "Hash = " << hashes << endl;
//            //if (hashes > 0)
//            //{
//            if ((double) L2.size() / (double) ReadSize > .6) {
//                ReadSize = L2.size();
//                lines++;
//                sequences.push_back(L2);
//                qual.push_back(L4);
//                string depths = "";
//
//                if (hashes > 0) {
//                    if (b[0] == 0) {
//                        strand.push_back(".");
//                    } else if (b[4] == 0) {
//                        strand.push_back("+");
//                    } else if (b[4] == 1) {
//                        strand.push_back("-");
//                    }
//                } else
//                    strand.push_back(".");
//
//                unsigned char C = 1;
//
//                for (int i = 0; i < L2.length(); i++) {
//                    depths += C;
//                }
//
//                depth.push_back(depths);
//
//            } else {
//                Rejects++;
//            }
//            //}
//        }
    }

    cout << endl;
    int NumReads = sequences.size();
    cout << "\nDone reading in \n		 Read in a total of " << NumReads + Rejects
         << " and rejected " << Rejects << endl;
    clock_t St, Et;
    float Dt;

    struct timeval start, end;
    gettimeofday(&start, NULL);
    int FoundMatch = 0;
    St = clock();

    for (std::vector<string>::size_type i = 0; i < sequences.size(); i++) {
        string A = sequences[i];
        string A_qualual = qual[i];
        string Adep = depth[i];
        string Astr = strand[i];

        if (FullOut) {
            cout
                    << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<**************************************************************>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
                    << endl;
        }

        Et = clock();

        if ((int) i % 10 < 1) {
            gettimeofday(&end, NULL);
            float Dt = end.tv_sec - start.tv_sec;
            cout << "aligning " << i << " of " << NumReads << ", \% done = "
                 << ((double) i / (double) NumReads) * 100.00 << ", TotalTime= " << Dt << " , second per read = "
                 << Dt / i << ", \% finding match = " << ((double) FoundMatch / (double) i) * 100.00 << "\r";
        }

        if (FullOut) {
            Dt = ((double) (Et - St)) / CLOCKS_PER_SEC;
            cout << "aligning " << i << " of " << NumReads
                 << "\% done = " << ((double) i / (double) NumReads) * 100.00
                 << ", TotalTime= " << Dt << " , second per read = " << Dt / i
                 << ", \% finding match = "
                 << ((double) FoundMatch / (double) i) * 100.00 << endl
                 << A << endl;

            for (int z = 0; z < Adep.length(); z++) {
                int bam = Adep.c_str()[z];
                cout << bam;
            }
            cout << endl;
        }

        int k = -1;
        int bestIndex = -1;
        bool PerfectMatch = false;

        int booya = Align3(sequences, qual, A, A_qualual, i, k, bestIndex, MinPercent, PerfectMatch, MinOverlap, Threads);
        if (FullOut) {
            cout << "best forward score is " << booya << " k is " << k << endl;
        }

        if (!(PerfectMatch)) {
            string revA = Util::RevComp(A);
            string revA_qualual = Util::RevQual(A_qualual);
            string revAdep = Util::RevQual(Adep);
            string revAstr = FlipStrands(Astr);
            int revk = -1;
            int revbestIndex = -1;
            int revbooya =
                    Align3(sequences, qual, revA, revA_qualual, i, revk, revbestIndex, MinPercent, PerfectMatch, MinOverlap,
                           Threads);
            if (FullOut) {
                cout << "best reverse score is " << revbooya << " k is " << revk
                     << endl;
            }

            if (revbooya > booya) {
                A = revA;
                A_qualual = revA_qualual;
                Adep = revAdep;
                Astr = revAstr;
                k = revk;
                booya = revbooya;
                bestIndex = revbestIndex;
            }

        } else {
            if (FullOut) {
                cout << "Perfect Match Found, Skipping Reverse Search" << endl;
            }
        }

        if (booya < MinOverlap) {
            if (FullOut) {
                cout << "No good match found, skipping" << endl;
            }
        } else {
            FoundMatch++;
            string B = sequences[bestIndex];
            string B_qualual = qual[bestIndex];
            string Bdep = depth[bestIndex];
            string Bstr = strand[bestIndex];

            if (k > 0) {

                if (FullOut) {
                    cout << "found match at " << k << endl;

                    for (int z = 0; z < k; z++) {
                        cout << "+";
                    }

                    cout << A << endl << B << endl;

                    for (int z = 0; z < Bdep.length(); z++) {
                        int bam = Bdep.c_str()[z];
                        cout << bam;
                    }

                    cout << endl;
                }
            } else {
                if (FullOut) {
                    cout << "found match at " << k << endl;
                    cout << A << endl;

                    for (int z = 0; z < abs(k); z++) {
                        cout << "-";
                    }

                    cout << B << endl;

                    for (int z = 0; z < abs(k); z++) {
                        cout << "-";
                    }

                    for (int z = 0; z < Bdep.length(); z++) {
                        int bam = Bdep.c_str()[z];
                        cout << bam;
                    }

                    cout << endl;
                }
            }

            if (i == bestIndex) {
                cout << "ERROR ____________________ SAME READS " << endl;
            }
            if (A.size() != Adep.size() && B.size() != Bdep.size()) {
                cout << " ERRPR somethis the wrong size\n	A= " << A.size()
                     << " Ad = " << Adep.size() << " B= " << B.size()
                     << " Bd = " << Bdep.size() << endl;
            }

            string combined =
                    ColapsContigs(A, B, k, A_qualual, B_qualual, Adep, Bdep, Astr, Bstr);

            if (combined.size() != Bdep.size()) {
                cout << " ERRPR combined is the wrong size\n	C= " << combined.size()
                     << " Bd = " << Bdep.size() << endl;
            }

            sequences[bestIndex] = combined;
            qual[bestIndex] = B_qualual;
            depth[bestIndex] = Bdep;
            strand[bestIndex] = Bstr;
            sequences[i] = "moved";

            if (FullOut) {
                cout << combined << endl;

                for (int z = 0; z < Bdep.length(); z++) {
                    int bam = Bdep.c_str()[z];
                    cout << bam;
                }
            }
        }
    }

    cout << "\n\nRESULTS\n";
    int count = 0;
    cout << "sequences size = " << sequences.size() << endl;
    for (int i = 0; i < sequences.size(); i++) {

        if (FullOut) {
            cout << ">>>>" << sequences[i] << endl;
        }

        if (sequences[i] != "moved" && sequences[i].size() >= 95) {
            string rDep = depth[i];
            int maxDep = -1;

            for (int z = 0; z < rDep.size(); z++) {
                unsigned char bam = rDep.c_str()[z];
                if ((int) bam > maxDep) {
                    maxDep = (int) bam;
                }
            }

            if (maxDep >= MinCoverage /*&& maxDep >= 2*/) {

                if (sequences[i].size() != qual[i].size() && qual[i].size() != depth[i].size()) {
                    cout << "ERROR, read "
                         << "@NODE_" << argv[6] << "_" << i << "_L=" << sequences[i].size()
                         << "_D=" << maxDep
                         << " Has the wrong size, Seq = " << sequences[i].size()
                         << " Qual = " << qual[i].size() << " Dep = " << depth[i].size()
                         << endl;
                }

                count++;
                int F = 0;
                int R = 0;
                compresStrand(strand[i], F, R);
                report << "@NODE_" << argv[6] << "_" << i << "_L=" << sequences[i].size() << "_D=" << maxDep << ":" << F
                       << ":" << R << ":" << endl;
                report << sequences[i] << endl;
                report << "+" << endl;
                report << qual[i] << endl;

                Depreport << "@NODE_" << argv[6] << "_" << i << "_L=" << sequences[i].size() << "_D=" << maxDep << ":"
                          << F << ":" << R << ":" << endl;
                Depreport << sequences[i] << endl;
                Depreport << "+" << endl;
                Depreport << qual[i] << endl;
                Depreport << strand[i] << endl;
                unsigned char C = depth[i].c_str()[0];
                int booya = C;
                Depreport << booya;

                for (int w = 1; w < depth[i].size(); w++) {
                    C = depth[i].c_str()[w];
                    booya = C;
                    Depreport << " " << booya;
                }

                Depreport << endl;
            }
        }
    }
    if (MinCoverage <= 1) {
        for (int i = 0; i < Unsequences.size(); i++) {

            if (FullOut) {
                cout << ">>>>" << Unsequences[i] << endl;
            }

            if (Unsequences[i] != "moved" && Unsequences[i].size() >= 95) {
                int maxDep = -1;

                if (Unsequences[i].size() != Unqual[i].size() &&
                    Unqual[i].size() != Undepth[i].size()) {
                    cout << "ERROR, read "
                         << "@NODE_" << argv[6] << "_" << i << "_L=" << Unsequences[i].size()
                         << "_D=" << maxDep
                         << " Has the wrong size, Seq = " << Unsequences[i].size()
                         << " Qual = " << Unqual[i].size() << " Dep = " << Undepth[i].size()
                         << endl;
                }

                count++;
                report << "@NODE_" << argv[6] << "_" << i << "_L=" << Unsequences[i].size()
                       << "_D" << maxDep << endl;
                report << Unsequences[i] << endl;
                report << "+" << endl;
                report << Unqual[i] << endl;

                Depreport << "@NODE_" << argv[6] << "_" << i
                          << "_L=" << Unsequences[i].size() << "_D" << maxDep << endl;
                Depreport << Unsequences[i] << endl;
                Depreport << "+" << endl;
                Depreport << Unqual[i] << endl;
                Depreport << Unstrand[i] << endl;
                unsigned char C = Undepth[i].c_str()[0];
                int booya = C;
                Depreport << booya;

                for (int w = 1; w < Undepth[i].size(); w++) {
                    C = Undepth[i].c_str()[w];
                    booya = C;
                    Depreport << " " << booya;
                }

                Depreport << endl;
            }
        }
    } else
        cout << "min coverage = " << MinCoverage << " skipping Unaligned sequences" << endl;
    cout << "\nWrote " << count << " sequences" << endl;
    report.close();
}
