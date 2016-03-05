#ifndef MAIN_H
#define MAIN_H

/**
 *  main.h
 *  MetaQuant
 *
 *  Kevin McLoughlin
 *  Based on code from eXpress, created by Adam Roberts in 2011.
 *  Copyright 2014 Kevin McLoughlin, Adam Roberts. All rights reserved.
 */

#include "config.h"
#include "logger.h"
#include <algorithm>
#include <limits>
#include <boost/foreach.hpp>
#include <cmath>
#include <cassert>
#include <limits>

#define foreach BOOST_FOREACH

class MapParser;
class TargetTable;
class BiasBoss;
class MismatchTable;
class LengthDistribution;

extern Logger logger;
/**
 * A global bool that is true when processing is still occuring.
 * This is primarily used to notify the bias update thread to stop running.
 */
extern bool running;
/**
 * A global bool that is true when the auxilary params are finished burning in.
 * This is primarily used to notify the bias update thread to stop updating
 * certain parameters.
 */
extern bool burned_out;
/**
 * A global bool that is true when edit detection is enabled
 */
extern bool edit_detect;
/**
 * A global size_t for the maximum allowed indel size.
 */
extern size_t max_indel_size;
/**
 *  Global defaults for alpha, beta priors in stick-breaking distribution
 */
extern double alpha_prior;
extern double beta_prior;
/**
 * An enum for the allowed directions of reads.
 *  BOTH - either direction is acceptable.
 *  FR - The first (or only read if single-end) must be mapped to the forward
 *       strand and the second to the reverse.
 *  RF - The first (or only read if single-end) must be mapped to the reverse
 *       strand and the second to the forward.
 *  F  - The single-end read must be mapped to the forward strand.
 *  R  - The single-end read must be mapped to the reverse strand.
 */
enum Direction { FR, RF, R, F, BOTH };
/**
 * A global variable specifying which direction(s) is (are) allowed for input
 * fragments.
 */
extern Direction direction;
/**
 * A global size_t specifying the maximum read length supported.
 */
extern size_t max_read_len;
/**
 * A global size_t specifying the number of possible nucleotides.
 */
const size_t NUM_NUCS = 4;
/**
 * A global character array specifying the nucleotide ordering and encoding.
 */
const char NUCS[] = {'A','C','G','T'};
/**
 * A global double representing the log of 0.
 */
const double LOG_0 = HUGE_VAL;
/**
 * A global double representing the log of 1.
 */
const double LOG_1 = 0;
/**
 * A global double representing the log of 0.25.
 */
const double LOG_QUARTER = log(0.25);
/**
 * A global double specifying the default epsilon value to be used in approx_eq.
 */
const double EPSILON = 0.000001;
const double LOG_EPSILON = log(EPSILON);
const double LOG_MAX = log(std::numeric_limits<double>::max());

/**
 * A global double specifying the minimum logged target abundance to be used
 * for updating the background bias model
 */
const double LOG_MIN_TARG_RHO = log(0.0001);
/**
 * Global function that determines if two doubles are within some epsilon of
 * each other.
 * @param a first double to be compared.
 * @param b second double to be compared.
 * @param eps the maximum allowed distance between the doubles for them to be
 *        considered approximately equal.
 * @return True iff a and b are within eps of each other.
 */
inline bool approx_eq(double a, double b, double eps=EPSILON) {
  return fabs(a-b) <= eps;
}

/**
 * Global function to calculate the log of the sum of 2 logged values
 * efficiently.
 * @param x a double for the first logged value in the sum.
 * @param y a double for the second logged value in the sum.
 * @return a double for the log of exp(x)+exp(y).
 */
inline double log_add(double x, double y) {
  if (fabs(x) == LOG_0) {
    return y;
  }
  if (fabs(y) == LOG_0) {
    return x;
  }

  if (y > x) {
    std::swap(x,y);
  }

  double sum = x+log(1+exp(y-x));
  return sum;
}
/**
 * Global function to calculate the log of the difference of 2 logged values
 * efficiently.
 * @param x a double for the logged minuend.
 * @param y a double for the logged subtrahend.
 * @return a double for the log of exp(x)-exp(y).
 */
inline double log_sub(double x, double y) {
  if (fabs(y) == LOG_0) {
    return x;
  }
  
  // Have to be careful of numerical issues, so we allow y to be slightly
  // greater than x.
  if (x <= y) {
    assert(approx_eq(x, y));
    return LOG_0;
  }
  

  double diff = x+log(1-exp(y-x));
  return diff;
}
/**
 * Global function to determine if a logged value is 0 in non-log space.
 * @param x a double for the logged value to be tested.
 * @return True if exp(x)==0 (without rounding).
 */
inline double islzero(double x)
{
  return (fabs(x) == LOG_0);
}

/**
 * Global function to exponentiate a logged value safely.
 * @param x a double for the logged value to be exponentiated.
 * @return exp(x) or 0 if x == log(0).
 */
inline double sexp(double x)
{
  if (islzero(x)) {
    return 0.0;
  }
  return exp(x);
}

#endif
