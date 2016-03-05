/**
 *  bundles.h
 *  MetaQuant
 *
 *  Kevin McLoughlin
 *  Based on code from eXpress, created by Adam Roberts in 2011.
 *  Copyright 2014 Kevin McLoughlin, Adam Roberts. All rights reserved.
 */

#ifndef metaquant_bundles_h
#define metaquant_bundles_h

#include <boost/thread.hpp>
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>
#include <vector>

class Target;
typedef size_t TargID;
typedef boost::unordered_map<size_t, float> CovarMap;

/**
 * The CovarTable is a sparse matrix for storing and updating pairwise
 * covariances between targets.
 **/
class CovarTable {
  /**
   * A private map to look up the covariance for pairs of Targets by their
   * combined hashed TargIDs. These values are stored positive and logged,
   * even though the true covariances are negative.
   */
  CovarMap _covar_map;

public:
  /**
   * CovarTable Constructor.
   */
  CovarTable() {};
  /**
   * A member function that increases the covariance between two targets by the
   * specified amount (logged). These values are stored positive even though
   * the true covariance is negative.
   * @param targ1 one of the targets in the pair.
   * @param targ2 the other target in the pair.
   * @param covar a double specifying the amount to increase the pair's
   *        covariance by (logged, positive).
   */
  void increment(TargID targ1, TargID targ2, double covar);
  /**
   * A member function that returns the covariance between two targets.
   * The returned value will be the the negative of the true value (logged).
   * @param targ1 one of the targets in the pair.
   * @param targ2 the other target in the pair.
   * @return The negative of the pair's covariance (logged).
   */
   double get(TargID targ1, TargID targ2);
  /**
   * A member function that returns the number of pairs of targets with non-zero
   * covariance.
   * @return The number of target pairs with non-zero covariance.
   */
  size_t size() const { return _covar_map.size(); }
};

class BundleTable;

/**
 * The Bundle class keeps track of a group of targets that have shared ambiguous
 * (multi-mapped) reads.
 **/
class Bundle {
  /**
   * A private vector that stores pointers to all targets in the bundle.
   */
  std::vector<Target*> _targets;
  /**
   * A private size_t that stores the total number of observed fragments mapped
   * to targets in the bundle.
   */
  size_t _counts;
  /**
   * A private double that stores the total mass of observed fragments mapped
   * to targets in the bundle (logged), including the initial pseudo-mass.
   */
  double _mass;
  /**
   * A private pointer to the bundle this one was merged into for threadsafe
   * collapsing. 
   */
  Bundle* _merged_into;
  /**
   * Mutex for threadafety.
   */
  mutable boost::mutex _mut;
  
  friend class BundleTable;

public:
  /**
   * Bundle Constructor.
   * @param targ a pointer to the initial Target object in the bundle.
   */
  Bundle(Target* targ);
  /**
   * A private method for returning the root of the merge tree that this bundle
   * is a node in.
   * @return A pointer to the bundle at the root of the merge tree for this 
   *         bundle.
   */
  const Bundle* get_rep() const;
  /**
   * A member function that increases the total bundle observed fragment counts
   * by a given amount.
   * @param incr_amt the amount to increase the counts by.
   */
  void incr_counts(size_t incr_amt=1);
  /**
   * A member function that increases the total bundle mass (logged)
   * by a given amount.
   * @param incr_amt the amount to increase the mass by (logged).
   */
  void incr_mass(double incr_amt);
  /**
   * A member function that resets the Bundle mass to (log) 0.
   * Call is not passed on to _merged_into.
   */
  void reset_mass();
  /**
   * An accessor for the number of Targets in the bundle.
   * @return The number of Targets in the bundle.
   */
  size_t size() const;
  /**
   * An accessor for a pointer to the vector of pointers to Targets in the
   * bundle. The returned value does not outlive this.
   * @return Pointer to the vector pointing to bundle Targets.
   */
  const std::vector<Target*>* targets() const { return &_targets; }
  /**
   * An accessor for the the total number of observed fragments mapped to
   * targets in the bundle.
   * @return The total number of fragments mapped to targets in the bundle.
   */
  size_t counts() const;
  /**
   * An accessor for the the total mass of observed fragments mapped to
   * targets in the bundle (logged), including the initial pseudo-mass.
   * @return The total mass of fragments mapped to targets in the bundle.
   */
  double mass() const;
};

typedef boost::unordered_set<Bundle*> BundleSet;

/**

 * The BundleTable class keeps track of the Bundle objects for a given run. It
 * has the ability to create, delete, and merge bundles.
 **/
class BundleTable {
  /**
   * A private unordered_set to store all of the bundles.
   */
  BundleSet _bundles;
  /**
   * A private boolean specifying if methods needs to be threadsafe.
   */
  bool _threadsafe_mode;
  /**
   * A private mutex for threadsafety.
   */
  mutable boost::mutex _mut;
  /**
   * A private method for returning the root of the merge tree that the given
   * bundle is a node in.
   * @return A pointer to the bundle at the root of the merge tree for the given
   *         bundle.
   */
  Bundle* get_rep(Bundle* b);
public:
  /**
   * BundleTable Constructor.
   */
  BundleTable();
  /**
   * BundleTable deconstructor.  Deletes all Bundle objects.
   */
  ~BundleTable();
  /**
   * A member function that returns the set of current Bundle objects. The
   * returned object does not outlive this.
   * @return A reference to the unordered_set containing all current Bundle
   *         objects.
   */
  const BundleSet& bundles() const { return _bundles; }
  /**
   * An accessor for the current number of Bundles.
   * @return The current number of Bundles.
   */
  size_t size() const { return _bundles.size(); }
  /**
   * A member function that creates a new Bundle, initially containing only the
   * single given Target.
   * @param targ a pointer to the only Target initially contained in the Bundle
   * @return A pointer to the new Bundle object
   */
  Bundle* create_bundle(Target* targ);
  /**
   * A member function that merges two Bundle objects into one. The Targets are
   * all moved to the larger bundles and the other is deleted.
   * @param b1 a pointer to one of the Bundle objects to merge.
   * @param b2 a pointer to the other Bundle object to merge.
   * @return A pointer to the merged Bundle object.
   */
  Bundle* merge(Bundle* b1, Bundle* b2);
  /**
   * Collapses the merge tree so that all targets are placed in the target list
   * of the root node and all other nodes are deleted.
   */
  void collapse();
  /**
   * Accessor for whether or not the BundleTable is in threadsafe mode.
   * @return True if the BundleTable is in threadsafe mode.
   */
  bool threadsafe_mode() const { return _threadsafe_mode; }
  /**
   * Mutator for threadsafe mode.
   * @param mode bool specifying if threadsafe mode should be enabled (true)
   *        or disabled (false)
   */
  void threadsafe_mode(bool mode) { _threadsafe_mode = mode; }
};

#endif
