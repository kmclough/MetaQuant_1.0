/**
 *  biascorrection.cpp
 *  MetaQuant
 *
 *  Kevin McLoughlin
 *  Based on code from eXpress, created by Adam Roberts in 2011.
 *  Copyright 2014 Kevin McLoughlin, Adam Roberts. All rights reserved.
 *
 */

#include <algorithm>
#include <boost/assign.hpp>
#include <cassert>

#include "biascorrection.h"

#include "fragments.h"
#include "frequencymatrix.h"
#include "main.h"
#include "sequence.h"
#include "targets.h"

using namespace std;

// size of bias window
const int WINDOW = 21;
// center position in window
const int CENTER = 11;
// number of bases on either side of center
const int SURROUND = 10;

SeqWeightTable::SeqWeightTable(size_t window_size, size_t order, double alpha)
    : _order(order),
      _observed(order, window_size, window_size, alpha),
      _expected(order, window_size, order+1, EPSILON) {
}

SeqWeightTable::SeqWeightTable(size_t window_size, size_t order,
                               string param_file_name, string identifier)
    : _order(order),
      _observed(order, window_size, window_size, 0),
      _expected(order, window_size, order+1, 0) {
        
  //TODO: Allow for orders and window sizes to be read from param file.
  
  ifstream infile (param_file_name.c_str());
  const size_t BUFF_SIZE = 99999;
  char line_buff[BUFF_SIZE];
  if (!infile.is_open()) {
    logger.severe("Unable to open parameter file '%s'.",
                  param_file_name.c_str());
  }

  do {
    infile.getline (line_buff, BUFF_SIZE, '\n');
  } while (strncmp(line_buff, identifier.c_str(), identifier.size()));
  infile.getline (line_buff, BUFF_SIZE, '\n');

  do {
    infile.getline (line_buff, BUFF_SIZE, '\n');
  } while (strcmp(line_buff, "\tObserved Conditional Probabilities"));
  infile.getline (line_buff, BUFF_SIZE, '\n');
  
  //FG
  for (size_t k = 0; k < (size_t)WINDOW; ++k) {
    infile.getline (line_buff, BUFF_SIZE, '\n');
    char *p = strtok(line_buff, "\t");
    for (size_t i=0; i < pow(4.0, (double)min(order, k)); ++i) {
      for (size_t j=0; j < NUM_NUCS; ++j) {
        p = strtok(NULL, "\t");
        _observed.update(k, i, j, log(strtod(p,NULL)));
      }
    }
  }
  
  infile.getline (line_buff, BUFF_SIZE, '\n');
  infile.getline (line_buff, BUFF_SIZE, '\n');
  //BG
  for (size_t k = 0; k < order+1; ++k) {
    infile.getline (line_buff, BUFF_SIZE, '\n');
    char *p = strtok(line_buff, "\t");
    for (size_t i=0; i < pow(4.0, (double)min(order, k)); ++i) {
      for (size_t j=0; j < NUM_NUCS; ++j) {
        p = strtok(NULL, "\t");
        _expected.update(k, i, j, log(strtod(p,NULL)));
      }
    }
  }
}


void SeqWeightTable::copy_observed(const SeqWeightTable& other) {
  _observed = other._observed;
}

void SeqWeightTable::copy_expected(const SeqWeightTable& other) {
  _expected = other._expected;
}

void SeqWeightTable::increment_expected(const Sequence& seq, double mass,
                                        const vector<double>& fl_cdf) {
  _expected.fast_learn(seq, mass, fl_cdf);
}

void SeqWeightTable::normalize_expected() {
  _expected.calc_marginals();
}

void SeqWeightTable::increment_observed(const Sequence& seq, size_t i,
                                        double mass) {
  int left = (int)i - SURROUND;
  _observed.update(seq, left, mass);
}

double SeqWeightTable::get_weight(const Sequence& seq, size_t i) const {
  int left = (int)i - SURROUND;
  return _observed.seq_prob(seq, left) - _expected.seq_prob(seq, left);
}

void SeqWeightTable::append_output(ofstream& outfile) const {
  char buff[200];
  string header = "";
  for (int i = 0; i < WINDOW; i++) {
    sprintf(buff, "\t%d", i-CENTER+1);
    header += buff;
  }
  header += '\n';

  outfile << "\tObserved Marginal Distribution\n" << header;

  for (size_t j = 0; j < NUM_NUCS; j++) {
    outfile << NUCS[j] << ":\t";
    for (int i = 0; i < WINDOW; i++) {
      outfile << scientific << sexp(_observed.marginal_prob(i,j)) << "\t";
    }
    outfile<<endl;
  }

  outfile << "\tObserved Conditional Probabilities\nPosition\t";

  for (size_t j = 0; j < pow((double)NUM_NUCS, (double)_order+1); j++) {
    string s = "->";
    s += NUCS[j & 3];
    size_t cond = j >> 2;
    for (size_t k = 0; k < _order; ++k) {
      s = NUCS[cond & 3] + s;
      cond = cond >> 2;
    }
    outfile << s << '\t';
  }
  outfile << endl;

  for (int i = 0; i < WINDOW; i++) {
    outfile << i-SURROUND << ":\t";
    for (size_t j = 0; j < pow((double)NUM_NUCS, (double)_order+1); j++) {
      outfile << scientific << sexp(_observed.transition_prob(i,j>>2,j&3))
              << "\t";
    }
    outfile<<endl;
  }

  outfile << "\tBackground Conditional Probabilities\nPosition\t";

  for (size_t j = 0; j < pow((double)NUM_NUCS, (double)_order+1); j++) {
    string s = "->";
    s += NUCS[j & 3];
    size_t cond = j >> 2;
    for (size_t k = 0; k < _order; ++k) {
      s = NUCS[cond & 3] + s;
      cond = cond >> 2;
    }
    outfile << s << '\t';
  }
  outfile << endl;

  for (size_t i = 0; i < _order+1; i++) {
    outfile << i+1 << ":\t";
    for (size_t j = 0; j <  pow((double)NUM_NUCS, (double)_order+1); j++) {
      outfile << scientific << sexp(_expected.transition_prob(i,j>>2,j&3))
              << "\t";
    }
    outfile<<endl;
  }
}

BiasBoss::BiasBoss(size_t order, double alpha)
    : _order(order),
      _5_seq_bias(WINDOW, order, alpha),
      _3_seq_bias(WINDOW, order, alpha){
}

BiasBoss::BiasBoss(size_t order, string param_file_name)
    : _order(order),
      _5_seq_bias(WINDOW, order, param_file_name, ">5"),
      _3_seq_bias(WINDOW, order, param_file_name, ">3") {
}

void BiasBoss::copy_observations(const BiasBoss& other) {
  _5_seq_bias.copy_observed(other._5_seq_bias);
  _3_seq_bias.copy_observed(other._3_seq_bias);
}

void BiasBoss::copy_expectations(const BiasBoss& other) {
  _5_seq_bias.copy_expected(other._5_seq_bias);
  _3_seq_bias.copy_expected(other._3_seq_bias);
}

void BiasBoss::update_expectations(const Target& targ, double mass,
                                   const vector<double>& fl_cdf) {
  if ((mass == LOG_0) || (mass < LOG_MIN_TARG_RHO)) {
    return;
  }

  if (direction != R) {
    const Sequence& seq_fwd = targ.seq(0);
    _5_seq_bias.increment_expected(seq_fwd, mass, fl_cdf);
  }
  if (direction != F) {
    const Sequence& seq_rev = targ.seq(1);
    _3_seq_bias.increment_expected(seq_rev, mass, fl_cdf);
  }
}

void BiasBoss::normalize_expectations() {
  _5_seq_bias.normalize_expected();
  _3_seq_bias.normalize_expected();
}

void BiasBoss::update_observed(const FragHit& hit, double normalized_mass)
{
  assert (hit.pair_status() != PAIRED || (int)hit.length() > WINDOW);

  const Sequence& t_seq_fwd = hit.target()->seq(0);
  const Sequence& t_seq_rev = hit.target()->seq(1);

  if (hit.pair_status() != RIGHT_ONLY) {
    _5_seq_bias.increment_observed(t_seq_fwd, hit.left(), normalized_mass);
  }
  if (hit.pair_status() != LEFT_ONLY) {
    _3_seq_bias.increment_observed(t_seq_rev, t_seq_rev.length()-hit.right(),
                                   normalized_mass);
  }
}

double BiasBoss::get_target_start_bias(const Target& targ, size_t pos) const {

  const Sequence& t_seq_fwd = targ.seq(0);
  return _5_seq_bias.get_weight(t_seq_fwd, pos);
}

double BiasBoss::get_target_end_bias(const Target& targ, size_t pos) const {

  const Sequence& t_seq_rev = targ.seq(1);
  return _3_seq_bias.get_weight(t_seq_rev, targ.length()-pos-1);
}

double BiasBoss::get_average_target_bias(std::map<const size_t, float>& start_bias,
                                         std::map<const size_t, float>& end_bias,
                                         const Target& targ) const {
  double tot_start = LOG_0;
  double tot_end = LOG_0;

  for (size_t i = 0; i < targ.length(); ++i) {
    double start_bias_i;
    double end_bias_i;
    try {
      start_bias_i = start_bias.at(i);
    } catch (out_of_range) {
      start_bias_i = get_target_start_bias(targ, i);
    }
    try {
      end_bias_i = end_bias.at(i);
    } catch (out_of_range) {
      end_bias_i = get_target_end_bias(targ, i);
    }
    tot_start = log_add(tot_start, start_bias_i);
    tot_end = log_add(tot_end, end_bias_i);
  }

  //double avg_start_bias = (tot_start ) - (log((double)targ.length()));
  //double avg_end_bias = (tot_end) - (log((double)targ.length()));
  double avg_bias = (tot_start + tot_end) - (2*log((double)targ.length()));
  //logger.info("Target %s avg biases: start = %f, end = %f, total = %f", targ.name().c_str(),
  //   avg_start_bias, avg_end_bias, avg_bias);
  assert(!isnan(avg_bias));
  return avg_bias;
}

void BiasBoss::update_target_bias(std::map<const size_t, float>& start_bias,
                                 std::map<const size_t, float>& end_bias,
                                 const Target& targ) const {

  const Sequence& t_seq_fwd = targ.seq(0);
  const Sequence& t_seq_rev = targ.seq(1);

  for (std::map<const size_t, float>::iterator it = start_bias.begin();
       it != start_bias.end(); ++it) {
    //=-=TODO: Shouldn't this be it->first???
    size_t pos = it->second;
    start_bias[pos] = _5_seq_bias.get_weight(t_seq_fwd, pos);
  }
  for (std::map<const size_t, float>::iterator it = end_bias.begin();
       it != end_bias.end(); ++it) {
    size_t pos = it->second;
    end_bias[targ.length()-pos-1] = _3_seq_bias.get_weight(t_seq_rev, pos);
  }
}

void BiasBoss::append_output(ofstream& outfile) const {
  outfile << ">5' Sequence-Specific Bias\n";
  _5_seq_bias.append_output(outfile);
  outfile << ">3' Sequence-Specific Bias\n";
  _3_seq_bias.append_output(outfile);
}
