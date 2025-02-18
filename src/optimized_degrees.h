#pragma once

#include <iostream>
#include <NTL/RR.h>
#include <string>
#include <vector>
#include <fstream>
#include "func.h"

using namespace std;
using namespace NTL;

std::vector<int> compute_min_multdepth_update(RR alpha, RR epsilon, int depth, long maxdeg, bool is_comp);
std::vector<int> compute_min_multdepth(RR alpha, RR epsilon, long maxdeg, bool is_comp);
