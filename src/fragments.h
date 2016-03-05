/**
 *  fragments.h
 *  MetaQuant
 *
 *  Kevin McLoughlin
 *  Based on code from eXpress, created by Adam Roberts in 2011.
 *  Copyright 2014 Kevin McLoughlin, Adam Roberts. All rights reserved.
 **/

#ifndef FRAGMENTS_H
#define FRAGMENTS_H

#include <string>
#include <vector>
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <cassert>
#include <api/BamAlignment.h>
#include "sequence.h"
#include "boost/scoped_ptr.hpp"

typedef size_t TargID;
struct Library;
class Target;
class TargetTable;
class Fragment;
class FragTarget;
typedef std::unordered_map<TargID, FragTarget*> FragTargetMap;
typedef std::pair<TargID, FragTarget*> FragTargetPair;

/**
 * PairStatus enum.
 * PAIRED denotes that the fragment mapping has both paired end reads.
 * LEFT_ONLY denotes that the single read is not reverse complemented => its
 *           left end is the left fragment end.
 * RIGHT_ONLY denotes that the single read is reverse complemented => its right
 *            end is the right fragment end.
 **/
enum PairStatus { PAIRED, LEFT_ONLY, RIGHT_ONLY };

/**
 * The Indel struct stores the information for a single insertion or deletion.
 **/
struct Indel {
  /**
   * A public size_t for the position of the Indel in the read. 0-based.
   */
  size_t pos;
  /**
   * A public size_t for the length of the Indel in the read.
   */
  size_t len;
  /**
   * Indel constructor
   */
  Indel(size_t p, size_t l) : pos(p), len(l) {}
};

/**
 * The ReadHit struct stores information for a single read alignment.
 */
struct ReadHit {
  /**
   * A public string for the SAM "Query Template Name" (fragment name)
   */
  std::string name;
  /**
   * A public bool specifying if this read was sequenced first according to the
   * SAM flag.
   */
  bool first;
  /**
   * A public bool specifying if this read was reverse complemented in its
   * alignment according to the SAM flag. This would also imply that the read
   * is the left end of the fragment.
   */
  bool reversed;
  /**
   * A public TargID for the target mapped to.
   */
  size_t targ_id;
  /**
   * A public size_t containing the 0-based leftmost coordinate mapped to in the
   * target.
   */
  size_t left;
  /**
   * A public size_t containing the position following the 0-based rightmost
   * coordinate mapped to in the target.
   */
  size_t right;
  /**
   * The read sequence.
   */
  SequenceFwd seq;
  /**
   * A public vector of Indel objects storing all insertions to the reference in
   * the read. Insertions are stored in read order.
   */
  std::vector<Indel> inserts;
  /**
   * A public vector of Indel objects storing all insertions to the reference in
   * the read. Deletions are stored in read order.
   */
  std::vector<Indel> deletes;
  /**
   * A public BamAlignment object storing the raw alignment information from
   * BamTools for the read. Only valid if BAM file is input.
   */
  BamTools::BamAlignment bam;
  /**
   * A public string storing the raw alignment information from for the read.
   * Only valid if SAM file is input.
   */
  std::string sam;
  /**
   * A public int containing the left position for the mate of the read.
   * -1 if single-end fragment. This is temporarily used to help find the mate
   * but is not used after.
   */
  int mate_l;
};

/**
 * The HitParams struct stores likelihood information for a single hit of a
 * fragment.
 */
struct HitParams {
  double align_likelihood;
  double full_likelihood;
  double posterior;
};

/**
 * The FragTarget class stores parameters related to all alignments from a
 * fragment to a given target.
 **/
class FragTarget {
public:
  /**
   * Public TargID for the target mapped to.
   */
  TargID targ_id;
  /**
   * Public index of the target in the target table.
   */
  size_t targ_index;
  /**
   * Public double representing the log of the posterior probability that the 
   * parent fragment arose from the associated target.
   */
  double posterior;
  /**
   * Public double holding the log of the current expectation of the target sampling
   * probability, up to a constant.
   */
  double expect_log_sample_prob;
  /**
   * Public double holding the log of the sum of the proabilities over all alignments.
   */
  double total_align_likelihood;

  /**
   * FragTarget default constructor.
   * @param tid TargID for a target hit by fragment.
   */
  FragTarget(TargID tid) : targ_id(tid), targ_index(0), 
     posterior(LOG_0), total_align_likelihood(LOG_0) { }
  /**
   * FragTarget constructor.
   * @param tid TargID for a target hit by fragment.
   * @param targ_index index of the target in the target table.
   */
  FragTarget(TargID tid, size_t targ_ind) : targ_id(tid), targ_index(targ_ind), 
     posterior(LOG_0), total_align_likelihood(LOG_0) { }
};

/**
 * The FragHit struct stores the information for a single fragment alignment.
 **/
class FragHit {
  /**
   * Private pointer to the target mapped to.
   */
  Target* _target;
  /**
   * Private pointer to data for the upstream (left) read alignment (if it
   * exists). Pointer is deleted with this.
   */
  boost::scoped_ptr<ReadHit> _read_l;
  /**
   * Private pointer to data for the downstream (right) read alignment (if it
   * exists). Pointer is deleted with this.
   */
  boost::scoped_ptr<ReadHit> _read_r;

  /**
   * Private structure holding parameters specific to a FragHit
   */
  HitParams _params;
  
public:
  /**
   * FragHit constructor for single-end read.
   * @param h pointer to the ReadHit struct for the single-end read.
   */
  FragHit(ReadHit* h) : _target(NULL) {
    if (h->reversed) {
      _read_r.reset(h);
    } else {
      _read_l.reset(h);
    }
  }
  /**
   * Fraghit constructor for paired-end read.
   * @param l pointer to the ReadHit struct for the upstream (left) read.
   * @param r pointer to the ReadHit struct for the downstream (right) read.
   */
  FragHit(ReadHit* l, ReadHit* r) : _target(NULL), _read_l(l), _read_r(r) {
    assert(!l->reversed);
    assert(r->reversed);
    assert(l->name == r->name);
    assert(l->targ_id == r->targ_id);
    assert(l->left <= r->left);
    assert(l->first != r->first);
  }
  /**
   * Accessor for the name of the fragment.
   * @return The name of the fragment.
   */
  std::string frag_name() const {
    if (_read_l) {
      return _read_l->name;
    }
    assert(_read_r);
    return _read_r->name;
  }
  /**
   * Accessor for the hit parameters (likelihood, etc.).
   * @return A pointer to the hit parameters.
   */
  HitParams* params() { return &_params; }
  /**
   * Const accessor for the hit parameters (likelihood, etc.).
   * @return A const pointer to the hit parameters.
   */
  const HitParams* params() const { return &_params; }
  /**
   * Accessor for a pointer to the Target object the fragment is aligned to.
   * @return A pointer to the Target aligned to.
   */
  Target* target() const {
    return _target;
  }
  /**
   * Mutator for a pointer to the Target object the fragment is aligned to.
   * @param target a pointer to the Target aligned to.
   */
  void target(Target* target) {
    _target = target;
  }
  /**
   * Function to calculate the portion of the log likelihood specific to a
   * fragment - target alignment.
   * @param frag A const reference to the Fragment involved in the alignment.
   * @param t A pointer to the Target involved in the alignment.
   * @return The alignment portion of the log likelihood.
   */
  double align_likelihood(const Fragment& frag, Target* t) const;
  /**
   * Accessor for the ID of the target the fragment is aligned to.
   * @return The ID of the target aligned to.
   */
  TargID target_id() const {
    if (_read_l) {
      return _read_l->targ_id;
    }
    assert(_read_r);
    return _read_r->targ_id;
  }
  /**
   * Accessor for the leftmost position aligned to (0-based).
   * @return The leftmost position aligned to in the target.
   */
  size_t left() const {
    if (_read_l) {
      return _read_l->left;
    }
    assert(_read_r);
    return _read_r->left;
  }
  /**
   * Accessor for one position past the rightmost position aligned to (0-based).
   * @return One past the rightmost position aligned to in the target.
   */
  size_t right() const {
    if (_read_r) {
      return _read_r->right;
    }
    assert(_read_l);
    return _read_l->right;
  }
  /**
   * Accessor for the length of the fragment alignment. Returns 0 if the
   * fragment is single-end.
   * @return Length of fragment mapping.
   */
  size_t length() const {
    if (_read_l && _read_r) {
      return _read_r->right - _read_l->left;
    }
    return 0;
  }
  /**
   * Const accessor for the alignment of the read at the leftmost (5') end of
   * the fragment in target coordinates or NULL if it was not sequenced.
   * @return A const pointer to the 5' read alignment.
   */
  const ReadHit* left_read() const {
    if (_read_l) {
      return _read_l.get();
    }
    return NULL;
  }
  /**
   * Const accessor for the alignment of the read at the rightmost (3') end of
   * the fragment in target coordinates or NULL if it was not sequenced.
   * @return A const pointer to the 3' read alignment.
   */
  const ReadHit* right_read() const {
    if (_read_r) {
      return _read_r.get();
    }
    return NULL;
  }
  /**
   * Const accessor for the alignment of the first (or only) read sequenced in
   * the fragment.
   * @return A const pointer to the alignment of the first read.
   */
  const ReadHit* first_read() const {
    if (_read_l && _read_l->first) {
      return _read_l.get();
    }
    assert(_read_r);
    return _read_r.get();
  }
  /**
   * Const accessor for the alignment of the second read sequenced in
   * the fragment. Returns NULL if single-end.
   * @return A const pointer to the alignment of the second read.
   */
  const ReadHit* second_read() const {
    if (_read_l && !_read_l->first) {
      return _read_l.get();
    } else if (_read_r && !_read_r->first) {
      return _read_r.get();
    } else {
      return NULL;
    }
  }
  /**
   * Accessor for the alignment of the read at the leftmost (5') end of
   * the fragment in target coordinates or NULL if it was not sequenced.
   * @return A pointer to the 5' read alignment.
   */
  ReadHit* left_read() {
    return const_cast<ReadHit*>(const_cast<const FragHit*>(this)->left_read());
  }
  /**
   * Accessor for the alignment of the read at the rightmost (3') end of
   * the fragment in target coordinates or NULL if it was not sequenced.
   * @return A pointer to the 3' read alignment.
   */
  ReadHit* right_read() {
    return const_cast<ReadHit*>(const_cast<const FragHit*>(this)->right_read());
  }
  /**
   * Accessor for the alignment of the first (or only) read sequenced in
   * the fragment.
   * @return A pointer to the alignment of the first read.
   */
  ReadHit* first_read() {
    return const_cast<ReadHit*>(const_cast<const FragHit*>(this)->first_read());
  }
  /**
   * Accessor for the alignment of the second read sequenced in
   * the fragment. Returns NULL if single-end.
   * @return A pointer to the alignment of the second read.
   */
  ReadHit* second_read() {
    return const_cast<ReadHit*>(const_cast<const FragHit*>(this)->second_read());
  }
  /**
   * A member function returning whether the mapping is PAIRED, LEFT_ONLY, or
   * RIGHT_ONLY, as defined in the PairStatus enum definition.
   * @return The pair status of the mapping.
   */
  PairStatus pair_status() const {
    if (_read_l && _read_r) {
      return PAIRED;
    }
    if (_read_l) {
      return LEFT_ONLY;
    }
    return RIGHT_ONLY;
  }
};

/**
 * The Fragment class stores information for all alignments of a single fragment.
 * By design, only paired-end mappings of paired-end reads will be accepted.
 * All mappings of single-end reads will be accepted.
 **/
class Fragment {
  /**
   * A private vector of FragHit pointers containing all mappings of the
   * fragment.
   */
  std::vector<FragHit*> _frag_hits;
  /**
   * A private map from TargIDs to FragTarget objects representing targets
   * hit by this Fragment.
   */
  FragTargetMap _frag_target_map;
  /**
   * A private vector of ReadHit pointers containing single read mappings whose
   * pairs have not been found. Temporarily useful when parsing the input to
   * find mates but is not used after.
   */
  std::vector<ReadHit*> _open_mates;
  /**
   * A private string for the SAM "Query Template Name" (fragment name).
   */
  std::string _name;
  /**
   * A private double for the mass of the Fragment as determined by the
   * forgetting factor during processing.
   */
  double _mass;
  /**
   * A private double for the learning rate associated with the Fragment, used
   * in the online variational Bayes algorithm.
   */
  double _learning_rate;
  /**
   * A private pointer to the global variables associated with the library
   * this fragment is from.
   */
  Library* _lib;
  /**
   * A private method that searches for the mate of the given read mapping.
   * If found, the mates are combined into a single FragHit and added to
   * frag_hits. If not found, the read mapping is added to open_mates.
   * @om pointer to single-end read mapping to find the mate of.
   */
  void add_open_mate(ReadHit* om);

public:
  /**
   * Fragment Constructor.
   */
  Fragment(Library* lib);
  /**
   * Fragment destructor deletes all FragHit, ReadHit and FragTarget objects pointed to by
   * the Fragment.
   */
  ~Fragment();
  /**
   * Accessor for the global variables associated with the library this fragment
   * is from. Pointer outlives this.
   */
  const Library* lib() const { return _lib; }
  /**
   * A member function that adds a new ReadHit to the Fragment. If it is the 
   * first ReadHit, it sets the Fragment name. If the fragment is not paired, a
   * FragHit is created and added to _frag_hits. Otherwise, add_open_mate is
   * called.
   * @param r a pointer to the ReadHit to be added.
   * @return True iff the read name matches the Fragment name or it is the first
   *         read.
   */
  bool add_map_end(ReadHit* r);
  /**
   * A member function that adds a FragHit pointer to the end of the _frag_hits
   * vector, creates an associated FragTarget object if none exists already
   * for the FragHit's TargID, and if so links it to the Fragment through its
   * FragHitMap.
   */
  void add_frag_hit(FragHit* h);
  /**
   * A member function that returns the number of FragTarget objects linked to this fragment
   */
  size_t num_frag_targets() const {
    return(_frag_target_map.size());
  }
  /**
   * A member function that returns a pointer to the FragTarget object for this
   * fragment and the given TargID
   */
  FragTarget* get_frag_target(TargID targ_id) const {
    return(_frag_target_map.at(targ_id));
  }
  /**
   * A member function that performs the variational Bayes parameter updates. First
   * assigns posterior probabilities for each target hit by this fragment based on
   * the alignment likelihoods and the current target sampling probability
   * estimates; then updates the stick-breaking beta distribution parameters for 
   * each target hit by this fragment, based on the posteriors.
   */
  void var_bayes_updates();
  /**
   * A member function that returns a reference to the "Query Template Name".
   * @return Reference to the SAM "Query Template Name" (fragment name).
   */
  const std::string& name() const { return _name; }
  /**
   * An accessor for the number of valid alignments of the fragment.
   * @return Number of valid alignments for fragment.
   */
  size_t num_hits() const { return _frag_hits.size(); }
  /**
   * An accessor for a pointer to the FragHit at the given index.
   * @param i index of the FragHit requested.
   * @return A pointer to the FragHit at the given index.
   */
  FragHit* operator[](size_t i) const { return _frag_hits[i]; }
  /**
   * Accessor for the FragHit objects associated with the fragment. Returned
   * value does not outlive this.
   * @return Reference to a vector containing pointers to the FragHits.
   */
  const std::vector<FragHit*>& hits() const { return _frag_hits; }
  /**
   * A member function that returns a single FragHit of the fragment sampled at
   * random based on the probabalistic assignments. Returned value does not
   * outlive this.
   * @return A randomly sampled FragHit.
   */
  const FragHit* sample_hit() const;
  /**
   * Mutator for the forgetting mass of the fragment.
   * @param m a double representing the value to set the forgetting mass to.
   */
  void mass(double m) { _mass = m; }
  /**
   * An accessor for the forgetting mass of the fragment.
   * @return The forgetting mass of the fragment.
   */
  double mass() const { return _mass; }
  /**
   * Mutator for the learning rate of the fragment.
   * @param lrate a double representing the value to set the learning rate to.
   */
  void learning_rate(double lrate) { _learning_rate = lrate; }
  /**
   * An accessor for the learning rate of the fragment according to the forgetting
   * factor.
   * @return The learning rate of the fragment.
   */
  double learning_rate() const { return _learning_rate; }
  /**
   * A member function that sorts the FragHits by the TargID of the targets they
   * are aligned to.
   */
  void sort_hits();
  /**
   * An accessor that returns true iff the Fragment has paired alignments.
   * @return True iff the Fragment has paired alignments.
   */
  bool paired() const {
    if (_frag_hits.empty()) {
      return false;
    }
    return (_frag_hits[0]->pair_status() == PAIRED);
  }
};

#endif
