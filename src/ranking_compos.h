#pragma once

#include "utils-basics.h"
#include "utils-eval.h"
#include "utils-matrices.h"
#include "utils-ptxt.h"
#include "utils-remez.h"
#include <cassert>
#include <omp.h>
#include "math/chebyshev.h"
#include "schemerns/rns-scheme.h"
#include "scheme/ckksrns/ckksrns-cryptoparameters.h"
#include <iostream>
#include <NTL/RR.h>
#include <cmath>
#include "optimized_degrees.h"
#include <vector>
#include <functional>
#include <cstdlib>  // For system()
#include <string>   // For std::string and std::to_string
#include <sstream>  // For std::ostringstream
#include <utility>  // For std::pair
#include <fstream>    // For std::ifstream

using namespace std;
using namespace NTL;




Ciphertext<DCRTPoly> rankComposCipher(
    Ciphertext<DCRTPoly> c,
    const std::vector<std::vector<double>> coefficient_list,
    const size_t vectorLength,
    double leftBoundC,
    double rightBoundC,
    std::vector<int> M_degs,
    float alpha,
    std::string alpha_string,
    bool cmpGt = false
);




std::vector<Ciphertext<DCRTPoly>> rankComposCipherVector(
    const std::vector<Ciphertext<DCRTPoly>> &c,
    const size_t subVectorLength,
    double leftBoundC,
    double rightBoundC,
    uint32_t degreeC,
    bool cmpGt,
    bool complOpt,
    std::vector<int> M_degs
);


// fractional rank
std::vector<double> rankComposPlain(
    const std::vector<double> &vec,
    double epsilon = 0.0
);
