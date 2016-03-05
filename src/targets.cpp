/*
 *  targets.cpp
 *  MetaQuant
 *
 *  Kevin McLoughlin
 *  Based on code from eXpress, created by Adam Roberts in 2011.
 *  Copyright 2014 Kevin McLoughlin, Adam Roberts. All rights reserved.
 */

#include "main.h"
#include "targets.h"
#include "lengthdistribution.h"
#include "fragments.h"
#include "biascorrection.h"
#include "mismatchmodel.h"
#include "mapparser.h"
#include "sequence.h"
#include "library.h"
#include <iostream>
#include <fstream>
#include <cassert>
#include <stdio.h>
#include <limits>
#include <float.h>
#include <boost/thread.hpp>
#include <boost/chrono.hpp>
#include <boost/atomic.hpp>
// =-= For memory debugging only
#include "memdebug.h"

using namespace std;

Target::Target(TargID id, const std::string& name, const std::string& description,
               const std::string& seq,
               bool prob_seq, double alpha, double beta, const Librarian* libs,
               const BiasBoss* known_bias_boss, const LengthDistribution* known_fld)
   : _libs(libs),
     _id(id),
     _name(name),
     _description(description),
     _seq_f(seq, 0, prob_seq),
     _seq_r(_seq_f),
     _alpha_prior(alpha),
     _beta_prior(beta),
     _ret_params(&_curr_params),
     _uniq_counts(0),
     _tot_counts(0),
     _avg_bias(0),
     _avg_bias_buffer(0),
     _solvable(false) {
  if ((_libs->curr_lib()).bias_table) {
    _start_bias.reset(new std::map<const size_t, float>);
    _start_bias_buffer.reset(new std::map<const size_t, float>);
    _end_bias.reset(new std::map<const size_t, float>);
    _end_bias_buffer.reset(new std::map<const size_t, float>);
  }
  update_target_bias_buffer(known_bias_boss, known_fld);
  swap_bias_parameters();
}

void Target::round_reset() {
  // =-= TODO: Note that this only makes sense for batch VBayes. For online VB _ret_params 
  // =-= should always point to _curr_params. Reimplement this when we implement batch
  // =-= VBayes. Will need to implement a copy constructor for RoundParams, to handle
  // =-= the atomic spin_lock variable.
  //_last_params = _curr_params;
  //_curr_params = RoundParams();
  //_ret_params = &_last_params;
}

double Target::stick_breaking_frac() const {
  return digamma(_ret_params->alpha) - digamma(_ret_params->alpha + _ret_params->beta);
}

double Target::stick_breaking_complement() const {
  return digamma(_ret_params->beta) - digamma(_ret_params->alpha + _ret_params->beta);
}

double Target::stick_breaking_frac_var() const {
  int ierr = 0;
  return trigamma(_ret_params->alpha, &ierr) - trigamma(_ret_params->alpha + _ret_params->beta, &ierr);
}

double Target::stick_breaking_complement_var() const {
  int ierr = 0;
  return trigamma(_ret_params->beta, &ierr) - trigamma(_ret_params->alpha + _ret_params->beta, &ierr);
}

void Target::update_alpha_beta(double frag_posterior, double posterior_sum,
       double learning_rate, size_t total_frags) {
  while (_curr_params.spin_lock.exchange(1, boost::memory_order_acquire) == 1) {
    // spin until free
  }
  _curr_params.alpha = _curr_params.alpha * (1. - learning_rate) 
                       + learning_rate*(_alpha_prior + total_frags*sexp(frag_posterior));
  _curr_params.beta = _curr_params.beta * (1. - learning_rate)
                + learning_rate*(_beta_prior + total_frags*sexp(posterior_sum));
  _curr_params.spin_lock.store(0, boost::memory_order_release);
}

double Target::rho() const {
  return _relative_abundance;
}

void Target::add_hit(double p, double m) {
  _curr_params.mass = log_add(_curr_params.mass, p+m);
  if (p != LOG_1) {
    if (p != LOG_0) {
      _curr_params.ambig_mass = log_add(_curr_params.ambig_mass, p+m);
      _curr_params.tot_ambig_mass = log_add(_curr_params.tot_ambig_mass, m);
    }
  }
  (_libs->curr_lib()).targ_table->update_total_fpb(m - _cached_eff_len);
}

double Target::est_effective_length(const LengthDistribution* fld,
                                    bool with_bias) const {
  if (!fld) {
    fld = (_libs->curr_lib()).fld.get();
  }

  double eff_len = LOG_0;

  double log_length = log((double)length());
  if (log_length < fld->mean()) {
    eff_len = log_length;
  } else {
    for(size_t l = fld->min_val(); l <= min(length(), fld->max_val()); l++) {
      eff_len = log_add(eff_len, fld->pmf(l)+log((double)length()-l+1));
    }
  }
  
  if (with_bias) {
    eff_len += _avg_bias;
  }

  return eff_len;
}

double Target::start_bias(const BiasBoss* bias_table, const size_t start_pos) {
  double bias;
  try {
    bias = _start_bias->at(start_pos);
  } catch (out_of_range err) {
    bias = bias_table->get_target_start_bias(*this, start_pos);
    (*_start_bias)[start_pos] = bias;
  }
  return bias;
}

double Target::end_bias(const BiasBoss* bias_table, const size_t end_pos) {
  double bias;
  try {
    bias = _end_bias->at(end_pos);
  } catch (out_of_range err) {
    bias = bias_table->get_target_end_bias(*this, end_pos);
    (*_end_bias)[end_pos] = bias;
  }
  return bias;
}

double Target::cached_effective_length(bool with_bias) const {
  if (with_bias) {
    return _cached_eff_len + _avg_bias;
  }
  return _cached_eff_len;
}

void Target::update_target_bias_buffer(const BiasBoss* bias_table,
                                       const LengthDistribution* fld) {
  long update_time = 0;
  long efflen_time = 0;

  if (bias_table) {
    boost::chrono::high_resolution_clock::time_point update_start = 
       boost::chrono::high_resolution_clock::now();

    bias_table->update_target_bias(*_start_bias_buffer, *_end_bias_buffer, *this);

    boost::chrono::high_resolution_clock::time_point update_end = 
       boost::chrono::high_resolution_clock::now();
    update_time = (long)((update_end-update_start).count());

    // =-= REMOVE: avg_bias never gets used in calculating effective length
    //_avg_bias_buffer = bias_table->get_average_target_bias(*_start_bias_buffer, 
    //   *_end_bias_buffer, *this);

  }
  //assert(!isnan(_avg_bias_buffer) && !isinf(_avg_bias_buffer));
  boost::chrono::high_resolution_clock::time_point efflen_start = 
       boost::chrono::high_resolution_clock::now();
  _cached_eff_len_buffer = est_effective_length(fld, false);
  boost::chrono::high_resolution_clock::time_point update_end = 
       boost::chrono::high_resolution_clock::now();

  efflen_time = (long)((update_end-efflen_start).count());

  //logger.info("Updated bias for target %s with hits at %d positions, rho = %f, update_time=%ld efflen_time=%ld",
  //    _name.c_str(), (*_start_bias_buffer).size(), sexp(_relative_abundance),
  //    update_time, efflen_time);
}

void Target::swap_bias_parameters() {
  _cached_eff_len = _cached_eff_len_buffer;
  _avg_bias = _avg_bias_buffer;
  _start_bias.swap(_start_bias_buffer);
  _end_bias.swap(_end_bias_buffer);
}


TargetTable::TargetTable(string targ_db_name, 
                         bool prob_seqs, bool known_aux_params, 
                         const Librarian* libs)
    :  _libs(libs) {
  string info_msg = "Loading target sequences";
  const Library& lib = _libs->curr_lib();
  const TargIndex& targ_index = lib.map_parser->targ_index();
  const TargIndex& targ_lengths = lib.map_parser->targ_lengths();
  const TargIDSet& hit_target_ids = lib.map_parser->hit_target_ids();
  if (lib.bias_table && !known_aux_params) {
    info_msg += " and measuring bias background";
  }
  info_msg += "...";
  logger.info(info_msg.c_str());

  boost::unordered_set<string> target_names;

  // Copy target index entries from the parser for targets with hits only into the
  // TargID set and vector objects
  size_t index = 0;
  for(TargIndex::const_iterator it = targ_index.begin(); it != targ_index.end(); ++it) {
    TargID targ_id = it->second;
    if (hit_target_ids.find(targ_id) != hit_target_ids.end()) {
      _hit_target_vec.push_back(targ_id);
      _targ_index_map[targ_id] = index++;
    }
  }

  // Find the length of the longest target sequence we might have to deal with
  size_t max_targ_length = 0;
  for(TargIndex::const_iterator it = targ_lengths.begin(); it != targ_lengths.end(); ++it) {
    size_t len = it->second;
    if (len > max_targ_length) {
      max_targ_length = len;
    }
  }
      
  // Instead of reading the whole FASTA file, we preprocess the reference FASTA
  // using code from refseq_db.py to generate a catalog, and then use fseek and fread
  // to only grab the sequences that have hits.

  string targ_catalog_file = targ_db_name + "_catalog.txt";
  string targ_data_file = targ_db_name + "_seq.dat";

  _target_db_fp = fopen(targ_data_file.c_str(), "rb");

  ifstream infile (targ_catalog_file.c_str());
  const size_t BUFF_SIZE = 2000;
  char line_buff[BUFF_SIZE];
  string name = "";
  long int offset = 0;
  size_t targ_len = 0;
  string targ_desc = "";
  size_t ntargs = 0;

  if (infile.is_open()) {
    char* seq_buff = (char*) malloc(max_targ_length*sizeof(char));
    while (infile.good()) {
      infile.getline(line_buff, BUFF_SIZE, '\n');
      // Split line into tab separated fields
      char *p = strtok(line_buff, "\t");
      int nf = 0;
      while (p && nf < 4) {
        switch(nf++) {
          case 0: {
            name = p;
            break;
          }
          case 1: {
            offset = (long int)atol(p);
            break;
          }
          case 2: {
            targ_len = (size_t)atoi(p);
            break;
          }
          case 3: {
            targ_desc = p;
            break;
          }
        }
        p = strtok(NULL, "\t");
      }
      long targ_size = add_targ(name, offset, targ_len, targ_desc, prob_seqs, known_aux_params, 
         alpha_prior, beta_prior, targ_index, targ_lengths, hit_target_ids,
         seq_buff);
      if (targ_size > 0) {
        ntargs++;
        //total_targ_data += targ_size;
        if (ntargs % 1000 == 1) {
          logger.info("Loaded target %d (%s) %s", ntargs, name.c_str(), 
              targ_desc.c_str());
        }
      }
    }
    infile.close();
    free(seq_buff);

    if (lib.bias_table && !known_aux_params) {
      lib.bias_table->normalize_expectations();
    }
  } else {
    logger.severe("Unable to open target DB catalog file '%s'.",
                  targ_catalog_file.c_str());
  }

  if (size() == 0) {
    logger.severe("No targets found in target DB catalog file '%s'.",
                  targ_catalog_file.c_str());
  }

  for (size_t i=0; i < size(); i++) {
    TargID targ_id = _hit_target_vec[i];
    if (_targ_map.find(targ_id) == _targ_map.end()) {
      logger.severe("Target %d not found in target DB catalog file "
                    "'%s'.", targ_id, targ_catalog_file.c_str());
    }
  }
  fclose(_target_db_fp);
  logger.info("Initialized %d targets.", ntargs);
  
}

TargetTable::~TargetTable() {
  foreach(TargPair targ_pair, _targ_map) {
    Target* targ = targ_pair.second;
    delete targ;
  }
}

long TargetTable::add_targ(const string& name, long int offset, size_t targ_len, 
                           const string& description, bool prob_seq,
                           bool known_aux_params, double alpha, double beta,
                           const TargIndex& targ_index,
                           const TargIndex& targ_lengths, 
                           const TargIDSet& hit_target_set,
                           char* seq_buffer) {

  long targ_size = 0;
  TargIndex::const_iterator it = targ_index.find(name);
  if (it == targ_index.end()) {
    logger.warn("Target '%s' exists in target DB but not alignment "
                   "(SAM/BAM) file.", name.c_str());
    return targ_size;
  }
  TargID targ_id = it->second;

  if (hit_target_set.find(targ_id) == hit_target_set.end()) {
    // Do nothing; we only add table entries for targets with hits
    return targ_size;
  }

  if (targ_lengths.find(name)->second != targ_len) {
    logger.severe("Target '%s' differs in length between target DB and "
                  "alignment (SAM/BAM) files (%d  vs. %d).", name.c_str(),
                  targ_len, targ_lengths.find(name)->second);
  }

  // Load the target sequence from the target DB data file
  int ferr = fseek(_target_db_fp, offset, SEEK_SET);
  if (ferr != 0) {
    logger.severe("Error trying to seek in target DB for target %s at offset %d",
       name.c_str(), offset);
  }
  size_t nread = fread((void *)seq_buffer, sizeof(char), targ_len, _target_db_fp);
  if (nread != targ_len) {
    logger.severe("Attempt to read sequence of length %d for target %s returned only %d chars",
       targ_len, name.c_str(), nread);
  }
  string seq(seq_buffer, targ_len);


  const Library& lib = _libs->curr_lib();
  const BiasBoss* known_bias_boss = (known_aux_params) ? lib.bias_table.get()
                                                       : NULL;
  const LengthDistribution* known_fld = (known_aux_params) ? lib.fld.get()
                                                           : NULL;
  
  Target* targ = new Target(targ_id, name, description, seq, prob_seq, alpha, beta, _libs,
                            known_bias_boss, known_fld);

  targ_size = sizeof(Target) + sizeof(SequenceFwd) + sizeof(SequenceRev) + 
     2*sizeof(RoundParams) + name.size() + description.size() + targ_len;
  if (lib.bias_table && !known_aux_params) {
    (lib.bias_table)->update_expectations(*targ);
  }
  _targ_map[targ_id] = targ;
  targ->bundle(_bundle_table.create_bundle(targ));

  return targ_size;
}

Target* TargetTable::get_targ(TargID id) {
    try {
      Target* t = _targ_map.at(id);
      return t;
    } catch (out_of_range err) {
      return NULL;
    }
}

Bundle* TargetTable::merge_bundles(Bundle* b1, Bundle* b2) {
  if (b1 != b2) {
    return _bundle_table.merge(b1, b2);
  }
  return b1;
}

void TargetTable::round_reset() {
  foreach(TargPair targ_pair, _targ_map) {
    Target* targ = targ_pair.second;
    targ->round_reset();
  }
}

void project_to_polytope(vector<Target*> bundle_targ,
                         vector<double>& targ_counts, double bundle_counts) {
  vector<bool> polytope_bound(bundle_targ.size(), false);
  while (true) {
    double unbound_counts = 0;
    double bound_counts = 0;
    for (size_t i = 0; i < bundle_targ.size(); ++i) {
      Target& targ = *bundle_targ[i];

      if (targ_counts[i] > targ.tot_counts()) {
        targ_counts[i] = targ.tot_counts();
        polytope_bound[i] = true;
      } else if (targ_counts[i] < targ.uniq_counts()) {
        targ_counts[i] = targ.uniq_counts();
        polytope_bound[i] = true;
      }

      if (polytope_bound[i]) {
        bound_counts += targ_counts[i];
      } else {
        unbound_counts += targ_counts[i];
      }
    }

    if (approx_eq(unbound_counts + bound_counts, bundle_counts)) {
      return;
    }
	
	  if (unbound_counts == 0) {
      polytope_bound = vector<bool>(bundle_targ.size(), false);
	    unbound_counts = bound_counts;
	    bound_counts = 0;
	  }

    double normalizer = (bundle_counts - bound_counts)/unbound_counts;
    for (size_t i = 0; i < bundle_targ.size(); ++i) {
      if (!polytope_bound[i]) {
        targ_counts[i] *= normalizer;
      }
    }
  }
}


void TargetTable::estimate_abundances() {
  double sum_rho = 0.0;
  for (size_t i = 0; i < size(); i++) {
    TargID targ_i = _hit_target_vec[i];
    Target* ti = get_targ(targ_i);
    double log_theta = ti->stick_breaking_frac();
    double var_log_theta = ti->stick_breaking_frac_var();
    for (size_t j = 0; j < i; j++) {
      TargID targ_j = _hit_target_vec[j];
      Target* tj = get_targ(targ_j);
      log_theta += tj->stick_breaking_complement();
      var_log_theta += tj->stick_breaking_complement_var();
    }
    ti->_sampling_prob = log_theta;
    // Estimate unnormalized log relative abundances and the standard error of the estimate
    ti->_relative_abundance = ti->_sampling_prob - ti->_cached_eff_len;
    // =-= Debugging code
    if (isnan(ti->_relative_abundance)) {
      logger.info("Target %s rho is NaN; sampling_prob = %f, cached_eff_len = %f, stick breaking frac = %f", ti->name().c_str(), ti->_sampling_prob, ti->_cached_eff_len, ti->stick_breaking_frac());
    }
    // =-= End debugging code

    sum_rho += sexp(ti->_relative_abundance);
    ti->_relative_abundance_sd = sqrt(var_log_theta);
  }
  // Normalize the relative abundances and confidence limits
  for (size_t i = 0; i < size(); i++) {
    TargID targ_i = _hit_target_vec[i];
    Target* ti = get_targ(targ_i);
    ti->_relative_abundance -= log(sum_rho);
    double log_ci = 1.959964 * ti->_relative_abundance_sd;
    ti->_relative_abundance_lo_cl = ti->_relative_abundance - log_ci;
    ti->_relative_abundance_hi_cl = ti->_relative_abundance + log_ci;
  }
}

void TargetTable::output_results(string output_dir, size_t tot_counts) {
  FILE * expr_file = fopen((output_dir + "/metaquant_results.txt").c_str(), "w");
  fprintf(expr_file, "bundle_id\ttarget_name\tlength\teff_length\t"
                     "rel_abundance\trel_abundance_conf_low\trel_abundance_conf_high\t"
                     "tot_counts\tuniq_counts\t"
                     "sampling_prob\tassigned_fragments\t"
                     "alpha\tbeta\tfpkm\tsolvable\tdescription\n");

  const double l_bil = log(1000000000.);
  
  size_t bundle_id = 0;
  
  bundle_id = 0;
  foreach (Bundle* bundle, _bundle_table.bundles()) {
    ++bundle_id;
    
    const vector<Target*>& bundle_targ = *(bundle->targets());
    
    for (size_t i = 0; i < bundle_targ.size(); ++i) {
      
      Target& targ = *bundle_targ[i];
      
      double fpkm = sexp(l_bil + targ._sampling_prob - targ._cached_eff_len);
      double eff_len = sexp(targ._cached_eff_len);
      double rho = sexp(targ._relative_abundance);
      double rho_lo_cl = sexp(targ._relative_abundance_lo_cl);
      double rho_hi_cl = sexp(targ._relative_abundance_hi_cl);
      double sample_prob = sexp(targ._sampling_prob);
      double assigned_frags = sample_prob * tot_counts;
      
      
      
      fprintf(expr_file, "" SIZE_T_FMT "\t%s\t" SIZE_T_FMT "\t%f\t"
              "%f\t%f\t%f\t"
              SIZE_T_FMT "\t" SIZE_T_FMT "\t" 
              "%f\t%.1f\t"
              "%f\t%f\t%f\t"
              "%c\t%s\n",
              bundle_id, targ.name().c_str(), targ.length(), eff_len,
              rho, rho_lo_cl, rho_hi_cl, 
              targ.tot_counts(), targ.uniq_counts(), 
              sample_prob, assigned_frags,
              targ.alpha(), targ.beta(), fpkm,
              (targ.solvable())?'T':'F', targ.description().c_str());
    }
  }
  fclose(expr_file);
}

double TargetTable::total_fpb() const {
  boost::unique_lock<boost::mutex>(_fpb_mut);
  return _total_fpb;
}

void TargetTable::update_total_fpb(double incr_amt) {
  boost::unique_lock<boost::mutex>(_fpb_mut);
  _total_fpb = log_add(_total_fpb, incr_amt);
}

void TargetTable::asynch_bias_update(boost::mutex* mutex, int sleep_time) {
  BiasBoss* bg_table = NULL;
  boost::scoped_ptr<BiasBoss> bias_table;
  boost::scoped_ptr<LengthDistribution> fld;

  bool burned_out_before = false;

  const Library& lib = _libs->curr_lib();

  while(running) {
    if (bg_table) {
      bg_table->normalize_expectations();
    }
    {
      boost::unique_lock<boost::mutex> lock(*mutex);
      logger.info("Synchronizing auxiliary parameter tables...");
      if(!fld) {
        fld.reset(new LengthDistribution(*(lib.fld)));
      } else {
        *fld = *(lib.fld);
      }
      if (lib.bias_table) {
        BiasBoss& lib_bias_table = *(lib.bias_table);
        if (!bias_table) {
          bias_table.reset(new BiasBoss(lib_bias_table));
        } else {
          lib_bias_table.copy_expectations(*bg_table);
          bg_table->copy_observations(lib_bias_table);
          bias_table.reset(bg_table);
        }
        bg_table = new BiasBoss(lib_bias_table.order(), 0);
      }
      // Update target abundance estimates. This is safe to do now since
      // process_fragment() is blocked by mutex for this code block.
      estimate_abundances();

      logger.info("Synchronized auxiliary parameter tables.");
    }

    if (!edit_detect && burned_out && burned_out_before) {
      break;
    }

    burned_out_before = burned_out;

    vector<double> fl_cdf = fld->cmf();

    // Buffer results of long computations.
    logger.info("Updating target bias buffers.");
    foreach(TargPair targ_pair, _targ_map) {
      Target* targ = targ_pair.second;
      // =-= Don't think we need to lock targets while calculating biases
      //targ->lock();
      targ->update_target_bias_buffer(bias_table.get(), fld.get());
      if (bg_table) {
        bg_table->update_expectations(*targ, targ->rho(), fl_cdf);
      }
      //targ->unlock();
    }
    logger.info("Updated target bias buffers.");
    {
      boost::unique_lock<boost::mutex> lock(*mutex);
      // Do quick atomic swap
      foreach(TargPair targ_pair, _targ_map) {
        Target* targ = targ_pair.second;
        targ->lock();
      }
      foreach(TargPair targ_pair, _targ_map) {
        Target* targ = targ_pair.second;
        targ->swap_bias_parameters();
        targ->unlock();
      }
    }
    logger.info("Swapped target bias parameters.");

    // Go to sleep for a while
    for (int remaining_sleep = sleep_time; remaining_sleep > 0; remaining_sleep -= 5) {
      if (!running) {
        break;
      }
      boost::this_thread::sleep_for(boost::chrono::seconds(5));
    }
  }

  if (bg_table) {
    delete bg_table;
  }
}
