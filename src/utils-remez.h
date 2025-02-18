#pragma once

#include "utils-basics.h"
#include "utils-eval.h"
#include "utils-matrices.h"
#include "ranking.h"
#include "minimum.h"
#include "sorting.h"
#include "openfhe.h"
#include "pke/cryptocontext.h"
#include "key/privatekey.h"
#include "key/publickey.h"
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

#include <cmath>
#include <chrono>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>


void run_remez(float approximation_param, std::vector<std::pair<RR, RR>> D, std::vector<int> degrees);

std::vector<double> read_coefficients(int degree, int i);

std::vector<std::vector<double>> read_coefficients_alpha(float alpha, std::vector<int> degrees);

float read_error();

float evaluateChebyshevBasis(const std::vector<double>& coefficients, float x);


float max_error_interval(const std::vector<double>& coefficients, float err);
float sign(float x);