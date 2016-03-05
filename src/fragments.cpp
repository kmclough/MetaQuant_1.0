/*
 *  fragments.cpp
 *  MetaQuant
 *
 *  Kevin McLoughlin
 *  Based on code from eXpress, created by Adam Roberts in 2011.
 *  Copyright 2014 Kevin McLoughlin, Adam Roberts. All rights reserved.
 */

#include "fragments.h"
#include "main.h"
#include "library.h"
#include "mismatchmodel.h"
#include "lengthdistribution.h"
#include "targets.h"
#include <string.h>
#include <stdlib.h>

using namespace std;

Fragment::Fragment(Library* lib) : _lib(lib) {}

Fragment::~Fragment() {
  for (size_t i = 0; i < num_hits(); i++) {
    delete _frag_hits[i];
  }

  // Delete FragTarget objects referenced by _frag_target_map
  foreach(FragTargetPair frag_targ_pair, _frag_target_map) {
    delete frag_targ_pair.second;
  }

  for (size_t i = 0; i < _open_mates.size(); i++) {
    delete _open_mates[i];
  }
}

bool Fragment::add_map_end(ReadHit* r)
{
  if (_name.empty()) {
    _name = r->name;
  } else if (_name != r->name) {
    //logger.info("add_map_end: name mismatch %s != %s", _name.c_str(), r->name.c_str());
    return false;
  }

  if (r->mate_l >= 0) {
    add_open_mate(r);
  } else {  // single-end fragment
    add_frag_hit(new FragHit(r));
  }

  return true;
}

void Fragment::add_frag_hit(FragHit* h) {
  _frag_hits.push_back(h);
  TargID targ_id = h->target_id();
  size_t targ_index = 0;
  if (!(_lib->targ_table.get() == 0)) {
    targ_index = _lib->targ_table->get_targ_index(targ_id);
  }
  FragTarget* frag_target = NULL;
  try {
    frag_target = _frag_target_map.at(targ_id);
  } catch (out_of_range) {
    frag_target = new FragTarget(targ_id, targ_index);
    _frag_target_map[targ_id] = frag_target;
  }
}

void Fragment::add_open_mate(ReadHit* nm) {
  bool found = false;

  for(vector<ReadHit*>::iterator it = _open_mates.begin();
      it != _open_mates.end(); ++it) {
    ReadHit* om = *it;
    if (nm->targ_id == om->targ_id &&
        (size_t)nm->mate_l == om->left &&
        (size_t)om->mate_l == nm->left &&
        nm->first != om->first &&
        nm->reversed != om->reversed) {
      FragHit* h = NULL;
      if (nm->left < om->left || (nm->left == om->left && om->reversed)) {
        h = new FragHit(nm, om);
      } else {
        h = new FragHit(om, nm);
      }

      found = true;
      add_frag_hit(h);
      _open_mates.erase(it);
      break;
    }
  }

  if (!found) {
    _open_mates.push_back(nm);
  }
}

bool fragtarget_compare(FragTarget* ft1, FragTarget* ft2) {
  return ft1->targ_index < ft2->targ_index;
}

void Fragment::var_bayes_updates() {

  // If fragment hits multiple targets, add target sampling probability to alignment 
  // likelihoods and normalize total likelihood across targets so likelihoods add to one.
  size_t num_targets = num_frag_targets();
  size_t i = 0;

  size_t ntargs = _lib->targ_table->size();

  if (num_targets > 1) {
    // Make an ordered list of targets hitting this fragment
    vector<FragTarget*> frag_targets(num_targets);
    i = 0;
    foreach(FragTargetPair frag_targ_pair, _frag_target_map) {
      FragTarget* frag_targ = frag_targ_pair.second;
      frag_targets[i++] = frag_targ;
    }
    sort(frag_targets.begin(), frag_targets.end(), fragtarget_compare);

    // Expectation step: Assign posterior probabilities to each target hit by fragment.
    double total_likelihood = LOG_0;
    // Compute expectation of log sampling probability for each target that hits fragment
    for (i = 0; i < num_targets; i++) {
      FragTarget* frag_targ = frag_targets[i];
      TargID targ_id = frag_targ->targ_id;
      Target* t = _lib->targ_table->get_targ(targ_id);

      frag_targ->expect_log_sample_prob = t->stick_breaking_frac();
      for (size_t j = 0; j < i; j++) {
        Target* u = _lib->targ_table->get_targ(frag_targets[j]->targ_id);
        frag_targ->expect_log_sample_prob += u->stick_breaking_complement();
      }
      frag_targ->posterior = frag_targ->total_align_likelihood + 
             frag_targ->expect_log_sample_prob;
      total_likelihood = log_add(total_likelihood, frag_targ->posterior);
    }
    // Scale the posteriors so they add to 1
    for (i = 0; i < num_targets; i++) {
      FragTarget* frag_targ = frag_targets[i];
      frag_targ->posterior -= total_likelihood;
    }

    // Maximization step: Update stick-breaking beta distribution parameters for each target
    for (size_t targ_ind = 0; targ_ind < ntargs; targ_ind++) {
      FragTarget* frag_targ = NULL;
      Target* targ_i = _lib->targ_table->get_targ_at_index(targ_ind);
      TargID targ_id = targ_i->id();
      double targ_i_posterior = LOG_0;
      try {
        frag_targ = _frag_target_map.at(targ_id);
        targ_i_posterior = frag_targ->posterior;
      } catch (out_of_range) {
      }
      double posterior_sum = LOG_0;
      for (size_t j = 0; j < num_targets; j++) {
        FragTarget* frag_targ_j = frag_targets[j];
        if (frag_targ_j->targ_index > targ_ind) {
          posterior_sum = log_add(posterior_sum, frag_targ_j->posterior);
        }
      }
      targ_i->update_alpha_beta(targ_i_posterior, posterior_sum, _learning_rate, 
        _lib->total_frags);
    }

  } else {
    // Expectation step

    // If fragment only hits one target, its posterior is 1
    FragTargetMap::iterator frag_targ_pair = _frag_target_map.begin();
    FragTarget* frag_targ = frag_targ_pair->second;
    frag_targ->posterior = LOG_1;
    size_t frag_targ_ind = frag_targ->targ_index;

    // Maximization step
    for (size_t targ_ind = 0; targ_ind < ntargs; targ_ind++) {
      Target* targ_i = _lib->targ_table->get_targ_at_index(targ_ind);
      double targ_i_posterior = LOG_0;
      double posterior_sum = LOG_0;
      if (targ_ind == frag_targ_ind) {
        targ_i_posterior = LOG_1;
      }
      if (frag_targ_ind > targ_ind) {
        posterior_sum = LOG_1;
      }
      targ_i->update_alpha_beta(targ_i_posterior, posterior_sum, _learning_rate, 
        _lib->total_frags);
    }
  }
}


const FragHit* Fragment::sample_hit() const {
  vector<double> probs(_frag_hits.size());
  probs[0] = sexp(_frag_hits[0]->params()->posterior);
  for (size_t i=1; i < _frag_hits.size(); ++i) {
    probs[i] = probs[i-1] + sexp(_frag_hits[i]->params()->posterior);
  }

  double r = rand()/double(RAND_MAX)*probs.back();
  size_t i = lower_bound(probs.begin(), probs.end(), r) - probs.begin();
  return _frag_hits[i];
}

bool fraghit_compare(FragHit* h1, FragHit* h2) {
  if (h1->target_id() == h2->target_id()) {
    return h1->left() < h2->left();
  } else {
    return h1->target_id() < h2->target_id();
  }
}

void Fragment::sort_hits() {
  sort(_frag_hits.begin(), _frag_hits.end(), fraghit_compare);
}

double FragHit::align_likelihood(const Fragment& frag, Target* t) const {

  const Library& lib = *frag.lib();

  double ll = LOG_1;

  const PairStatus ps = pair_status();

  if (lib.mismatch_table) {
    ll += (lib.mismatch_table)->log_likelihood(*this);
  }

  if (lib.bias_table) {
    if (ps != RIGHT_ONLY) {
      ll += t->start_bias(lib.bias_table.get(), left());
    }
    if (ps != LEFT_ONLY) {
      ll += t->end_bias(lib.bias_table.get(), right() - 1);
    }
  }
  
  if (ps == PAIRED) {
    ll += (lib.fld)->pmf(t->length());
  } else if (ps == LEFT_ONLY && t->length() - left() < (lib.fld)->max_val()) {
    ll += (lib.fld)->cmf(t->length() - left());
  } else if (ps == RIGHT_ONLY && right() < (lib.fld)->max_val()) {
    ll += (lib.fld)->cmf(right());
  }

  assert(!(isnan(ll)||isinf(ll)));
  return ll;
}

