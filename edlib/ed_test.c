//
// Created by Sanhu Li on 10/19/17.
//

#include <stdio.h>
#include <string.h>
#include "edlib.h"
#include "../utils.h"


int main() {
//    EdlibAlignResult result = edlibAlign("hello", 5, "world!", 6, edlibNewAlignConfig(42, EDLIB_MODE_HW, EDLIB_TASK_PATH,
//                                                                                      NULL, 0));
    const char *s1 = "AACTCCGACTGCTCAACAAGAGGTCACACCGTTAGGTCGACCTCAGCCCCGTACTGGCCGAAAGCGTGAGATGACACGGGGCAAGCTTGATGTTCCGAAC",
            *s2 = "AACTGCGACTGCTCAACAAGATGTCACACCGTTAGGTCGACCTCAGCCCCGTACTGGCCGAAAGCGTGAGATGACACGGGGCAAGCTTGATGTTCCGAAC";
    EdlibAlignResult result = edlibAlign(s1, strlen(s1), s2, strlen(s2), edlibNewAlignConfig(5000000, EDLIB_MODE_NW, EDLIB_TASK_PATH,
                                                                                            NULL, 0));
    if (result.status == EDLIB_STATUS_OK) {
    	printf("%d\n%s\n", result.editDistance, edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_STANDARD));
//        cout << result.editDistance << endl;
//        cout << edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_STANDARD) << endl;
    }
//    cout << "telephone" << endl;
//    cout << "-elephant" << endl;
    edlibFreeAlignResult(result);
//
//    std::vector<const char *> lines;
//    read_lines("S288C_reference_sequence_R64-2-1_20150113.fsa.chn", lines);
//    for (auto l : lines) {
//        cout << l << endl;
//    }
}