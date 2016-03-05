/**
 *  library.h
 *  MetaQuant
 *
 *  Kevin McLoughlin
 *  Based on code from eXpress, created by Adam Roberts in 2012.
 *  Copyright 2014 Kevin McLoughlin, Adam Roberts. All rights reserved.
 */

#ifndef metaquant_library_h
#define metaquant_library_h

#include <vector>
#include "boost/shared_ptr.hpp"

/**
 * The Library struct holds pointers to the global parameter tables for a set of
 * reads from the same library preparation.
 */
struct Library {
  /**
   * Path to the input file. Empty if streamed.
   */
  std::string in_file_name;
  /**
   * Path to the out file. Empty if alignments are not to be output.
   */
  std::string out_file_name;
  /**
   * A pointer to the MapParser for parsing the input alignment file for this
   * library.
   */
  boost::shared_ptr<MapParser> map_parser;
  /**
   * A pointer to the fragment length distribution object for this library.
   */
  boost::shared_ptr<LengthDistribution> fld;
  /**
   * A pointer to the MismatchTable containing the learned error distribution
   * for this library.
   */
  boost::shared_ptr<MismatchTable> mismatch_table;
  /**
   * A pointer to the BiasBoss containing the learned bias distribution for this
   * library. (optional)
   */
  boost::shared_ptr<BiasBoss> bias_table;
  /**
   * A pointer to the TargetTable containing the target parameters (abundance,
   * effective length) for this library.
   */
  boost::shared_ptr<TargetTable> targ_table;
  /**
   * The number of the next read to be processed (starting at 1).
   */
  size_t n;
  /**
   * The total number of fragments represented in the input alignment file.
   */
  size_t total_frags;
  /**
   * The mass of the next read to be processed (logged).
   */
  double mass_n;
  /**
   * Library constructor sets initial values for parameters
   */
  Library() : n(1), mass_n(0) {};
};

/**
 * The Librarian class keeps track of the different library objects for a run.
 */
class Librarian
{
  /**
   * A private vector of library structs for the different inputs in the run.
   */
  std::vector<Library> _libs;
  /**
   * The index of the library currently being processed.
   */
  size_t _curr;

public:
  /**
   * Librarian Constructor.
   * @param num_libs a size_t for the number of libraries to be processed in the
   *        run.
   */
  Librarian(size_t num_libs): _libs(num_libs), _curr(0) {}
  /**
   * An accessor for the Library struct at a given index. Returned value does
   * not outlive this.
   * @param i a size_t indexing the requested Library struct.
   * @return The Library struct at the given index.
   */
  Library& operator[](size_t i) {
        assert(i < _libs.size());
        return _libs[i];
  }
  /**
   * An accessor for the Library struct associated with the library currently
   * being processed. Returned value does not outlive this.
   * @return The Library struct indexed by _curr.
   */
  const Library& curr_lib() const { return _libs[_curr]; }
  /**
   * A mutator of the index of the library currently being processed.
   * @param i a size_t to set the index of the current Library struct to.
   */
  void set_curr(size_t i) {
    assert (i < _libs.size());
    _curr = i;
  }
  /**
   * An accessor for the number of Library structs. This should be equal to the
   * number of libraries to be processed in the run.
   * @return The number of Library structs.
   */
  size_t size() const { return _libs.size(); }
};

#endif
