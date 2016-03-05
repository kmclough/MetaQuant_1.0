/**
 *  markovmodel.cpp
 *  MetaQuant
 *
 *  Kevin McLoughlin
 *  Based on code from eXpress, created by Adam Roberts in 2012.
 *  Copyright 2014 Kevin McLoughlin, Adam Roberts. All rights reserved.
 */

#include "markovmodel.h"
#include "frequencymatrix.h"
#include "sequence.h"
#include "main.h"

using namespace std;

MarkovModel::MarkovModel(size_t order, size_t window_size,
                         size_t num_pos, double alpha)
    : _order((int)order),
      _window_size((int)window_size),
      _num_pos((int)num_pos),
      _params(num_pos, FrequencyMatrix<double>((size_t)pow((double)NUM_NUCS,
                                                           (double)order),
                                               NUM_NUCS, alpha)),
      _bitclear((1<<(2*order))-1) {
}

size_t MarkovModel::get_indices(const Sequence& seq, int left, vector<char>& indices) {
  assert(_order <= 3);
  int i = 0;
  int j = left;
  int seq_len = (int)seq.length();
  
  size_t cond = 0;
  
  if (left != 0 && left < _order) {
    i = _order-left;
    for (j=0; j < _order; j++) {
      cond = (cond << 2) + seq[j];
    }
  }
  
  indices = vector<char>(_window_size - i, -1);
  size_t start_index = i;
  
  while (i < _window_size && j < seq_len) {
    size_t index = min(i, _num_pos-1) - start_index;
    size_t curr = seq[j];
    cond = (cond << 2) + curr;
    indices[index] = (char)cond;
    cond &= _bitclear;
    i++;
    j++;
  }
  return start_index;
}

void MarkovModel::update(const Sequence& seq, int left, double mass) {
  int i = 0;
  int j = left;
  int seq_len = (int)seq.length();

  size_t cond = 0;
  
  if (left != 0 && left < _order) {
    i = _order-left;
    for (j=0; j < _order; j++) {
      cond = (cond << 2) + seq[j];
    }
  }

  while (i < _window_size && j < seq_len) {
    size_t index = min(i, _num_pos-1);
    size_t curr = seq[j];
    //if (seq.prob()) {
    //  for (size_t nuc = 0; nuc < NUM_NUCS; nuc++) {
    //    _params[index].increment(cond, nuc, seq.get_prob(j, nuc) + mass);
    //  }
    //} else {
      _params[index].increment(cond, curr, mass);
    //}
    cond = (cond << 2) + curr;
    cond &= _bitclear;
    i++;
    j++;
  }
}

void MarkovModel::update(size_t p, size_t i, size_t j, double mass) {
  _params[p].increment(i, j, mass);
}

vector<char> MarkovModel::get_indices(const Sequence& seq) {
  vector<char> indices(seq.length() - _order, -1);
  
  if (seq.length() < (size_t)_order) {
    return indices;
  }
  
  size_t cond = 0;
  for (int i = 0; i < _order; ++i) {
    cond = (cond << 2) + seq[i];
  }
  
  for (size_t i = _order; i < seq.length(); ++i) {
    size_t curr = seq[i];
    cond = (cond << 2) + curr;
    indices[i-_order] = (char)cond;
    cond &= _bitclear;
  }
  return indices;
}

void MarkovModel::fast_learn(const Sequence& seq, double mass,
                             const vector<double>& fl_cmf) {
  assert(_num_pos==_order+1);
  if (seq.length() < (size_t)_order) {
    return;
  }

  size_t cond = 0;
  for (int i = 0; i < _order; ++i) {
     cond = (cond << 2) + seq[i];
  }

  for (size_t i = _order; i < seq.length(); ++i) {
    size_t curr = seq[i];

    double mass_i = mass;
    if (seq.length()-i < fl_cmf.size()) {
      mass_i += fl_cmf[seq.length()-i];
    }

    if (seq.prob()) {
      for (size_t nuc = 0; nuc < NUM_NUCS; nuc++) {
        _params[_order].increment(cond, nuc, seq.get_prob(i, nuc) + mass_i);
      }
    } else {
      _params[_order].increment(cond, curr, mass_i);
    }
    cond = (cond << 2) + curr;
    cond &= _bitclear;
  }
}

void MarkovModel::calc_marginals() {
  assert(_num_pos==_order+1);
  for (int i = 0; i < _order; ++i) {
    for (size_t cond = 0; cond < pow((double)NUM_NUCS,
                                     (double)(_order)); cond++) {
      for (size_t nuc = 0; nuc < NUM_NUCS; ++nuc) {
        _params[i].increment(cond & ((1 << (2*i))-1) , nuc,
                             _params[_order]((cond << 2) + nuc, false));
      }
    }
  }
}

double MarkovModel::transition_prob(size_t p, size_t cond, size_t curr) const {
  assert(p < _params.size());
  return _params[p](cond, curr);
}

double MarkovModel::seq_prob(const Sequence& seq, int left) const {
  int i = 0;
  int j = left;
  int seq_len = (int)seq.length();

  size_t cond = 0;
  double v = 0;

  if (left < _order) {
    i = _order-left;
    for (j=0; j < min(_order, i); j++) {
      cond = (cond << 2) + seq[j];
    }
    v = i*LOG_QUARTER;
  }

  while (i < _window_size && j < seq_len) {
    size_t index = min(i, _num_pos -1);
    size_t curr = seq[j];
    if (seq.prob()) {
      double prob = LOG_0;
      for (size_t nuc = 0; nuc < NUM_NUCS; nuc++) {
        assert(!isnan(seq.get_prob(j, nuc)));
        assert(!isnan(_params[index](cond, nuc)));
        prob = log_add(prob, seq.get_prob(j, nuc) + _params[index](cond, nuc));
      }
      v += prob;
    } else {
      assert(!isnan(_params[index](cond, curr)));
      v += _params[index](cond, curr);
    }
    cond = (cond << 2) + curr;
    cond &= _bitclear;

    i++;
    j++;
  }
  assert(!isnan(v));
  return v;
}

double MarkovModel::marginal_prob(size_t w, size_t nuc) const {
  assert(w < _params.size());
  double marg = LOG_0;
  double tot = LOG_0;
  for (size_t cond = 0; cond < pow((double)NUM_NUCS,
                                   (double)(_order)); cond++) {
    marg = log_add(marg, _params[w]((cond << 2) + nuc, false));
    tot = log_add(tot, _params[w].sum(cond));
  }
  return marg-tot;
}
