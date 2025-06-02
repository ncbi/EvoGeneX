#ifndef __EVOGENEX_COMMON_H__
#define __EVOGENEX_COMMON_H__

#include <RcppEigen.h>
#include <iostream>
#include <limits>
#include <Eigen/Dense>

const double X_TOL = sqrt(std::numeric_limits<double>::epsilon());

#define KEEP_LOG 0

using namespace Rcpp;
using namespace std;
using namespace Eigen;

#endif // __EVOGENEX_COMMON_H__
