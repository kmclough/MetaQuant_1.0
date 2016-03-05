/**
 *  targets.h
 *  MetaQuant
 *
 *  Kevin McLoughlin
 *  Based on code from eXpress, created by Adam Roberts in 2011.
 *  Copyright 2014 Kevin McLoughlin, Adam Roberts. All rights reserved.
 */

#ifndef TARGETS_H
#define TARGETS_H

#include <boost/scoped_ptr.hpp>
#include "boost/shared_ptr.hpp"
#include <boost/thread.hpp>
#include <boost/unordered_map.hpp>
#include <boost/math/special_functions/digamma.hpp>
using boost::math::digamma;
#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <cstdio>
#include <vector>
#include "main.h"
#include "bundles.h"
#include "sequence.h"
#include "trigamma.h"
#include <boost/atomic.hpp>

class LengthDistribution;
class FragHit;
class BiasBoss;
class MismatchTable;
class Librarian;
class TargetTable;

/**
 * The RoundParams struct stores the target parameters unique to a given round
 * (iteration) of VBayes.
 **/
struct RoundParams {
  /**
   * A public double that stores the alpha parameter for the Beta distribution
   * governing the stick-breaking probability for this target
   */
  double alpha;
  /**
   * A public double that stores the beta parameter for the Beta distribution
   * governing the stick-breaking probability for this target
   */
  double beta;
  /**
   * A public double that stores the (logged) assigned mass based on posterior
   * fragment mapping probabilities.
   */
  double mass;
  /**
   * A public double that stores the (logged) assigned ambiguous mass based on
   * posterior fragment mapping probabilities.
   */
  double ambig_mass;
  /**
   * A public double that stores the (logged) total mass of ambiguous fragments
   * mapping to the target.
   */
  double tot_ambig_mass;
  /**
   * An atomic int variable that controls access to the alpha and beta variables
   * in this object.
   */
  boost::atomic<int> spin_lock;
  /**
   * RoundParams constructor sets initial values for parameters
   */
  RoundParams() : alpha(1.0), beta(1.0), mass(LOG_0), ambig_mass(LOG_0),
   tot_ambig_mass(LOG_0), spin_lock(0) {}
};

typedef size_t TargID;

/**
 * The Target class is used to store objects for the targets being mapped to.
 * Besides storing basic information about the object (id, length), it also
 * stores a mass based on the number of fragments mapping to the object as well
 * as parameters for variance. To help with updating these values, it computes
 * the likelihood that a given fragment originated from it. These values are
 * stored and returned in log space.
 **/
class Target {
  friend class TargetTable;
  /**
   * A private pointer to the struct containing pointers to the global
   * parameter tables (bias_table, mismatch_table, fld).
   */
  const Librarian* _libs;
  /**
   * A private TargID that stores the target's index in the TargetTable.
   */
  TargID _id;
  /**
   * A private string that stores the target name.
   */
  std::string _name;
  /**
   * A private string that stores the target description.
   */
  std::string _description;
  /**
   * A private Sequence object that stores the forward target sequence.
   */
  SequenceFwd _seq_f;
  /**
   * A private Sequence object that stores the reverse target sequence.
   */
  SequenceRev _seq_r;
  /**
   * A private double object that stores the alpha parameter for the prior
   * stick-breaking Beta distribution.
   */
  double _alpha_prior;
  /**
   * A private double object that stores the beta parameter for the prior
   * stick-breaking Beta distribution.
   */
  double _beta_prior;
  /**
   * A private RoundParams struct that stores the parameters for the current
   * online or batch VBayes round.
   */
  RoundParams _curr_params;
  /**
   * A private RoundParams struct that stores the parameters for the previous
   * batch VBayes round.
   */
  RoundParams _last_params;
  /**
   * A private pointer to the RoundParams that should be used in any
   * accessors.
   */
  RoundParams* _ret_params;
  /**
   * A private size_t that stores the number of fragments (non-logged)
   * uniquely mapping to this target.
   */
  size_t _uniq_counts;
  /**
   * A private size_t that stores the fragment counts (non-logged) for the
   * bundle. The total bundle counts is the sum of this value for all targets
   * in the bundle.
   */
  size_t _tot_counts;
  /**
   * A private pointer to the Bundle this Target is a member of.
   */
  Bundle* _bundle;
  /**
   * A private mutex to provide thread-safety for variables with threaded
   * update.
   */
  mutable boost::mutex _mutex;
  /**
   * A scoped pointer to a private float map storing the (logged) 5' bias
   * at each position with a read hit.
   */
  boost::scoped_ptr<std::map<const size_t, float> > _start_bias;
  /**
   * Buffers the start bias to allow for atomic updating.
   */
  boost::scoped_ptr<std::map<const size_t, float> > _start_bias_buffer;
  /**
   * A scoped pointer to a private float map storing the (logged) 3' bias
   * at each position with a read hit.
   */
  boost::scoped_ptr<std::map<const size_t, float> > _end_bias;
  /**
   * Buffers the end bias to allow for atomic updating.
   */
  boost::scoped_ptr<std::map<const size_t, float> > _end_bias_buffer;
  /**
   * A private double storing the (logged) product of the average 3' and 5'
   * biases for the target.
   */
  double _avg_bias;
  /**
   * Buffers the average bias to allow for atomic updating.
   */
  double _avg_bias_buffer;
  /**
   * A private double storing the most recent estimate of the log probability of
   * sampling a fragment from this target. Only updated at the end of each round.
   */
  double _sampling_prob;
  /**
   * A private double storing the most recent estimate of the log relative abundance
   * of this target. Only updated at the end of each round.
   */
  double _relative_abundance;
  /**
   * A private double storing the standard error of the most recent estimate of the log 
   * relative abundance of this target. Only updated at the end of each round.
   */
  double _relative_abundance_sd;
  /**
   * A private double storing the lower 95% confidence limit estimate of the log relative 
   * abundance of this target. Only updated at the end of each round.
   */
  double _relative_abundance_lo_cl;
  /**
   * A private double storing the upper 95% confidence limit estimate of the log relative 
   * abundance of this target. Only updated at the end of each round.
   */
  double _relative_abundance_hi_cl;
  /**
   * A private double storing the most recently updated (logged) effective
   * length as calculated by the bias updater thread.
   */
  double _cached_eff_len;
  /**
   * Buffers the cached effective length to allow for atomic updating.
   */
  double _cached_eff_len_buffer;
  /**
   * A private boolean specifying whether a unique solution exists. True iff
   * a unique read is mapped to the target or all other targets in a mapping
   * are solvable.
   */
  bool _solvable;

public:
  /**
   * Target Constructor.
   * @param id a unique TargID identifier.
   * @param name a string that stores the target name.
   * @param description a string that stores the target description.
   * @param seq a string that stores the target sequence.
   * @param prob_seq a bool that specifies if the sequence is to be treated
   *        probablistically. For RDD detection.
   * @param alpha a double that specifies the initial (prior) value for the alpha
   *        parameter in the target's stick-breaking Beta distribution.
   * @param beta a double that specifies the initial (prior) value for the beta
   *        parameter in the target's stick-breaking Beta distribution.
   * @param libs a pointer to the struct containing pointers to the global
   *        parameter tables (bias_table, mismatch_table, fld).
   * @param known_bias_boss a pointer to bias parameters provided as input, NULL
   *        if none given.
   * @param known_fld a pointer to a fragment length distribution provided as
   *        input, NULL if none given.
   */
  Target(TargID id, const std::string& name, const std::string& description, 
         const std::string& seq,
         bool prob_seq, double alpha, double beta, const Librarian* libs,
         const BiasBoss* known_bias_boss, const LengthDistribution* known_fld);
  /**
   * A member function that locks the target mutex to provide thread safety.
   * The lock should be held by any thread that calls a method of the Target.
   */
  void lock() const { _mutex.lock(); }
  /**
   * A member function that unlocks the target mutex.
   */
  void unlock() const { _mutex.unlock(); }
  /**
   * An accessor for the target name.
   * @return string containing target name.
   */
  const std::string& name() const { return _name; }
  /**
   * An accessor for the target description.
   * @return string containing target description.
   */
  const std::string& description() const { return _description; }
  /**
   * An accessor for the target id.
   * @return The target ID.
   */
  TargID id() const { return _id; }
  /**
   * An accessor for the the target's Sequence (const).
   * @param rev a bool specifying whether to return the reverse complement.
   * @return Const reference to the target's Sequence object.
   */
  const Sequence& seq(bool rev=false) const {
    if (rev) {
      return _seq_r;
    }
    return _seq_f;
  }
  /**
   * An accessor for the the target's Sequence (non-const).
   * @param rev a bool specifying whether to return the reverse complement.
   * @return Non-const reference to the target's Sequence object.
   */
  Sequence& seq(bool rev) {
    if (rev) {
      return _seq_r;
    }
    return _seq_f;
  }
  /**
   * An accessor for the length of the target sequence.
   * @return The target sequence length.
   */
  size_t length() const { return _seq_f.length(); }
  /**
   * An accessor for the current estimated log relative abundance of the target.
   * Only updated at the end of each round.
   * @return The current estimated rho.
   */
  double rho() const;
  /**
   * A member function that estimates the log of the stick-breaking 
   * fraction for this target, using the current alpha and beta parameters.
   */
  double stick_breaking_frac() const;
  /**
   * A member function that estimates the log of the complement of the stick-breaking 
   * fraction for this target, using the current alpha and beta parameters.
   */
  double stick_breaking_complement() const;
  /**
   * A member function that estimates the variance of the log of the stick-breaking 
   * fraction for this target, using the current alpha and beta parameters.
   */
  double stick_breaking_frac_var() const;
  /**
   * A member function that estimates the variance of the log of the complement of the 
   * stick-breaking fraction for this target, using the current alpha and beta parameters.
   */
  double stick_breaking_complement_var() const;
  /**
   * A member function that updates the alpha and beta parameters for the stick-breaking
   * Beta distribution of this target.
   * @param frag_posterior A double representing the log of the posterior likelihood
   *      for some fragment that hits this target.
   * @param posterior_sum A double representing the log of the sum of the posterior likelihoods
   *      for all higher numbered targets for the same fragment 
   * @param learning_rate A double representing the learning rate to be used in the
   *      parameter update.
   * @param total_frags A size_t representing the total number of fragments in the library.
   */
  void update_alpha_beta(double frag_posterior, double posterior_sum,
       double learning_rate, size_t total_frags);
  /**
   * A member function that increases the estimated fragment hit counts of the target
   * based on the posterior assignment probability of some FragHit.
   * @param p the posterior probability for the FragHit that is being added.
   * @param mass a double specifying the (logged) forgetting mass of the fragment being
   *        mapped.
   */
  void add_hit(double p, double mass);
  /**
   * An accessor for the current alpha parameter for the stick-breaking
   * Beta distribution of this target.
   * @return The alpha parameter.
   */
  double alpha() const { return _ret_params->alpha; }
  /**
   * An accessor for the current beta parameter for the stick-breaking
   * Beta distribution of this target.
   * @return The beta parameter.
   */
  double beta() const { return _ret_params->beta; }
  /**
   * An accessor for the current (logged) mass of fragments assigned to the target.
   * @return The logged mass.
   */
  double mass() const { return _ret_params->mass; }
  /**
   * An accessor for the the (logged) total mass of ambiguous fragments mapping
   * to the target.
   * @return The (logged) total mass of ambiguous fragments mapping to the
   *         target.
   */
  double tot_ambig_mass() const { return _ret_params->tot_ambig_mass; }
  /**
   * A member function that prepares the target object for the next round of
   * batch VBayes.
   */
  void round_reset();
  /**
   * An accessor for the current count of fragments mapped to this target
   * either uniquely or ambiguously.
   * @return The total fragment count.
   */
  size_t tot_counts() const { return _tot_counts; }
  /**
   * An accessor for the the current count of fragments uniquely mapped to this
   * target.
   * @return The unique fragment count.
   */
  size_t uniq_counts() const { return _uniq_counts; }
  /**
   * An accessor for the pointer to the Bundle this Target is a member of.
   * @return A pointer to the Bundle this target is a member of.
   */
  Bundle* bundle() const { return _bundle; }
  /**
   * A mutator to set the Bundle this Target is a member of.
   * @param b a pointer to the Bundle to set this Target as a member of.
   */
  void bundle(Bundle* b) { _bundle = b; }
  /**
   * A member function that increases the count of fragments mapped to this
   * target.
   * @param uniq a bool specifying whether or not the fragment uniquely maps
   *        to this target.
   * @param incr_amt a size_t to increase the counts by.
   */
  void incr_counts(bool uniq, size_t incr_amt = 1) {
    if (uniq) {
      _solvable = true;
    }
    _tot_counts += incr_amt;
    _uniq_counts += incr_amt * uniq;
  }
  /**
   * A member function that returns the log likelihood term due to the positional bias
   * for the start of a fragment. Returns a cached value if the target has one for the
   * given position; otherwise computes the value and caches it in the target's
   * _start_bias map.
   * @param bias_table Pointer to the bias table to be used for computing the bias, if necessary.
   * @param start_pos The position of the start of the fragment on the target sequence
   * @return The start bias portion of the log likelihood for a fragment-target alignment.
   */
  double start_bias(const BiasBoss* bias_table, const size_t start_pos) ;
  /**
   * A member function that returns the log likelihood term due to the positional bias
   * for the end of a fragment. Returns a cached value if the target has one for the
   * given position; otherwise computes the value and caches it in the target's
   * _end_bias map.
   * @param bias_table Pointer to the bias table to be used for computing the bias, if necessary.
   * @param end_pos The position of the end of the fragment on the target sequence
   * @return The end bias portion of the log likelihood for a fragment-target alignment.
   */
  double end_bias(const BiasBoss* bias_table, const size_t end_pos);
  /**
   * A member function that calculates and returns the estimated effective
   * length of the target (logged) using the average bias.
   * @param fld an optional pointer to a different LengthDistribution than the
   *       global one, for thread-safety.
   * @param with_bias a boolean specifying whether or not the average bias
   *        should be included in the return value.
   * @return The estimated effective length of the target calculated as
   *         \f$ \tilde{l} = \bar{bias}\sum_{l=1}^{L(t)} D(l)(L(t) - l + 1) \f$.
   */
  double est_effective_length(const LengthDistribution* fld = NULL,
                              bool with_bias=true) const;
  /**
   * An accessor for the most recently estimated effective length (logged) as
   * calculated by the bias updater thread.
   * @param with_bias a boolean specifying whether or not the average bias
   * should be included in the return value.
   * @return The cached effective length of the target.
   */
  double cached_effective_length(bool with_bias=true) const;
  /**
   * A member function that causes the target bias to be re-calculated by the
   * _bias_table based on curent parameters. The results are buffered until
   * swap_bias_parameters is called to allow for atomic updating.
   * @param bias_table a pointer to a BiasBoss to use as parameters. Bias not
   *        updated if NULL.
   * @param fld an optional pointer to a different LengthDistribution than the
   *        global one, for thread-safety.
   */
  void update_target_bias_buffer(const BiasBoss* bias_table = NULL,
                                 const LengthDistribution* fld = NULL);
  /**
   * Swaps in the buffered bias parameters for atomic updating. The target
   * mutex should be held by the caller.
   */
  void swap_bias_parameters();
  /**
   * An accessor for the _solvable flag.
   * @return a boolean specifying whether or not the target has a unique
   *         solution for its abundance estimate.
   */

  bool solvable() const { return _solvable; }
  /**
   * A mutator that sets the _solvable flag.
   * @param a boolean specifying whether or not the target has a unique solution
   *        for its abundance estimate.
   */
  void solvable(bool s) { _solvable = s; }
};


typedef boost::unordered_map<TargID, Target*> TargMap;
typedef std::pair<TargID, Target*> TargPair;
typedef boost::unordered_map<TargID, size_t> TargIndexMap;
typedef std::set<TargID> TargIDSet;
typedef std::vector<TargID> TargIDVector;
typedef boost::unordered_map<std::string, size_t> TargIndex;
typedef boost::unordered_map<size_t, float> CovarMap;

/**
 * The TargetTable class is used to keep track of the Target objects for a run.
 * The constructor parses a catalog file to generate the Target objects and stores
 * them in a map keyed by their string id.
 **/
class TargetTable {
  /**
   * A private pointer to the struct containing pointers to the global parameter
   * tables (bias_table, mismatch_table, fld).
   */
  const Librarian* _libs;
  /**
   * A private file pointer for the target database data file
   */
  FILE* _target_db_fp;
  /**
   * A private map to look up pointers to Target objects by their TargID id.
   */
  TargMap _targ_map;
  /**
   * A private vector for iterating over targets with at least one fragment hit
   */
  TargIDVector _hit_target_vec;
  /*
   * A private map from TargIDs to indices in _hit_target_vec
   */
  TargIndexMap _targ_index_map;
  /**
   * The private table to keep track of Bundle objects.
   */
  BundleTable _bundle_table;
  /**
   * A private table to look up the covariance for pairs of Targets by their
   * combined hashed TargIDs. These values are stored logged and positive, even
   * though they are negative.
   */
  CovarTable _covar_table;
  /**
   * A private double that stores the (logged) total mass per base
   * (including pseudo-counts) to allow for rho calculations.
   * TODO: Do we still need this?
   */
  double _total_fpb;
  /**
   * A private mutex to make accesses to _total_fpb thread-safe.
   */
  mutable boost::mutex _fpb_mut;

  /**
   * A private function that validates and adds a target pointer to the table.
   * @param name the name of the trancript.
   * @param offset the offset of the target sequence in the target DB data file.
   * @param targ_len the length of the target sequence.
   * @param description a string describing the target sequence.
   * @param prob_seqs a bool that specifies if the sequence is to be treated
   *        probablistically, for RDD detection.
   * @param known_aux_params a bool that is true iff the auxiliary parameters
   *        (fld, bias) are provided and need not be learned.
   * @param alpha a double that specifies the prior alpha for the stick-breaking
   *        Beta distribution of the target (non-logged).
   * @param beta a double that specifies the prior beta for the stick-breaking
   *        Beta distribution of the target (non-logged).
   * @param targ_index the target-to-index map from the alignment file.
   * @param targ_lengths the target-to-length map from the alignment file, for
   *        validation.
   * @param hit_target_set a set of TargIDs for targets with hits in the alignment file.
   * @return The size of the target sequence data loaded for this target, if it was loaded,
   *         or zero if the target was not loaded because it had no hits.
   */
  long add_targ(const std::string& name, long int offset, size_t targ_len, 
                const std::string& description, bool prob_seqs,
                bool known_aux_params, double alpha, double beta,
                const TargIndex& targ_index, const TargIndex& targ_lengths,
                const TargIDSet& hit_target_set, char* seq_buffer);

public:
  /**
   * TargetTable Constructor.
   * @param targ_db_name a string storing the path and filename prefix for the
   *        target database files
   * @param prob_seqs a bool that specifies if the sequence is to be treated
   *        probablistically, for RDD detection.
   * @param known_aux_params a bool that is true iff the auxiliary parameters
   *        (fld, bias) are provided and need not be learned.
   * @param libs a pointer to the struct containing pointers to the global
   *        parameter tables (bias_table, mismatch_table, fld).
   */
  TargetTable(std::string targ_db_name, 
              bool prob_seqs, bool known_aux_params, 
              const Librarian* libs);
  /**
   * TargetTable Destructor. Deletes all of the target objects in the table.
   */
  ~TargetTable();
  /**
   * A member function that returns a pointer to the target with the given id.
   * @param id of the target queried.
   * @return A pointer to the target with the given id.
   */
  Target* get_targ(TargID id);
  /**
   * A member function that returns a pointer to the vector of TargIDs for all targets with hits
   */
  TargIDVector* get_hit_target_ids() { return &_hit_target_vec; }
  /**
   * A member function that returns the index of the given TargID in the vector of targets
   * with hits.
   * @param targ_id The TargID to be located in the hit target vector.
   * @return The index of the given TargID in the hit target vector.
   */
  size_t get_targ_index(TargID targ_id) { return _targ_index_map.at(targ_id); }
  /**
   * A member function that readies all Target objects in the table for the next
   * round of batch VBayes.
   */
  void round_reset();
  /**
   * An accessor for the number of targets with hits in the table.
   * @return The number of targets with hits in the table.
   */
  size_t size() const { return _hit_target_vec.size(); }
  /**
   * An accessor for the Target at a given index
   * @param i the index of the target requested.
   * @return A pointer to the Target at the requested index.
   */
  Target* get_targ_at_index(const size_t i) {
    assert(i < _hit_target_vec.size());
    return _targ_map[_hit_target_vec[i]];
  }
  /**
   * An accessor for the (logged) total mass per base, including pseudo-counts.
   * @return The (logged) total mass per base, including pseudo-counts.
   */
  double total_fpb() const;
  /**
   * a member function that increments the (logged) total mass per base.
   * @param incr_amt the (logged) amount to increment by.
   */

  void update_total_fpb(double incr_amt);
  /**
   * A member function that increases the (logged) covariance between two
   * targets by the specified amount. These values are stored positive even
   * though they are negative.
   * @param targ1 one of the targets in the pair
   * @param targ2 the other target in the pair
   * @param covar a double specifying the amount to increase the pair's
   *        covariance by (logged)
   */
  void update_covar(TargID targ1, TargID targ2, double covar) {
    _covar_table.increment(targ1, targ2, covar);
  }
  /**
   * An accessor for the covariance between two targets. These returned value
   * will be the log of the negative of the true value.
   * @param targ1 one of the targets in the pair.
   * @param targ2 the other target in the pair.
   * @return The negative of the pair's covariance (logged).
   */
  double get_covar(TargID targ1, TargID targ2) {
    return _covar_table.get(targ1, targ2);
  }
  /**
   * An accessor for number of pairs of targets with non-zero covariance.
   * @return The number of target pairs with non-zero covariance.
   */
  size_t covar_size() const { return _covar_table.size(); }
  /**
   * A member function that merges the given Bundles.
   * @param b1 a pointer to the first Bundle to merge.
   * @param b2 a pointer to the second Bundle to merge.
   * @return A pointer to the merged Bundle.
   */
  Bundle* merge_bundles(Bundle* b1, Bundle* b2);
  /**
   * An accessor for the number of bundles in the partition.
   * @return The number of bundles in the partition.
   */
  size_t num_bundles() const { return _bundle_table.size(); }
  /**
   * Calculate normalized relative abundances for each target from stick-breaking parameters
   * and effective lengths.
   */
  void estimate_abundances();
  /**
   * A member function that outputs the final expression data in a file called
   * 'metaquant_results.txt' in the given output directory.
   * @param output_dir the directory to output the expression file to.
   * @param tot_counts the total number of observed mapped fragments.
   */
  void output_results(std::string output_dir, size_t tot_counts);
  /**
   * A member function to be run asynchronously that continuously updates the
   * background bias values, target bias values, and target effective lengths.
   * @param mutex a pointer to the mutex to be used to protect the global fld
   *        and bias tables during updates.
   */
  void asynch_bias_update(boost::mutex* mutex, int sleep_time);
  void enable_bundle_threadsafety() { _bundle_table.threadsafe_mode(true); }
  void disable_bundle_threadsafety() { _bundle_table.threadsafe_mode(false); }
  /**
   * Collapses the merge trees in the BundleTable.
   */
  void collapse_bundles() { _bundle_table.collapse(); }
};

#endif
