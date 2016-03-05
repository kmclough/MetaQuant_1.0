/**
 *  mapparser.h
 *  MetaQuant
 *
 *  Kevin McLoughlin
 *  Based on code from eXpress, created by Adam Roberts in 2011.
 *  Copyright 2014 Kevin McLoughlin, Adam Roberts. All rights reserved.
 */

#ifndef metaquant_mapparser_h
#define metaquant_mapparser_h

#include <api/BamReader.h>
#include <api/BamWriter.h>
#include <boost/scoped_ptr.hpp>
#include <boost/thread.hpp>
#include <boost/unordered_map.hpp>
#include <string>
#include <vector>

#include <iostream>

class Fragment;
class TargetTable;
class FragHit;
struct ParseThreadSafety;
struct ReadHit;
struct Library;

typedef size_t TargID;
typedef boost::unordered_map<std::string, TargID> TargIndex;
typedef std::set<TargID> TargIDSet;

/**
 * The Parser class is an abstract class for implementing a SAMParser or
 * BAMParser. It fills Fragment objects by parsing an input file in SAM/BAM
 * format.
 **/
class Parser {
 protected:
  /**
   * The private target-to-index map.
   */
  TargIndex _targ_index;
  /**
   * The private target-to-length map.
   */
  TargIndex _targ_lengths;
  /**
   * A private set of TargIDs for targets with hits from the reads being parsed.
   */
  TargIDSet _hit_target_ids;
  /**
   * A private pointer to the current/last read mapping being parsed.
   */
  ReadHit* _read_buff;

 public:
  /**
   * Dummy destructor.
   */
  virtual ~Parser(){};
  /**
   * An accessor for the SAM header string.
   * @return The SAM header string.
   */
  virtual const std::string header() const=0;
  /**
   * An accessor for the target name to index map. Returns a reference that
   * does not outlive this.
   * @return Reference to the target-to-index map.
   */
  const TargIndex& targ_index() const { return _targ_index; }
  /**
   * An accessor for the target-to-length map. Returns a reference that does
   * not outlive this.
   * @return Reference to the target-to-length map.
   */
  const TargIndex& targ_lengths() const { return _targ_lengths; }
  /**
   * An accessor for the set of target IDs with hits. Returns a reference that does
   * not outlive this.
   * @return Reference to the hit target ID set.
   */
  TargIDSet& hit_target_ids() { return _hit_target_ids; }
  /**
   * A member function that loads all mappings of the next fragment into the
   * given Fragment object.
   * @param f the empty Fragment to add mappings to.
   * @return True iff more reads remain in the SAM/BAM file/stream.
   */
  virtual bool next_fragment(Fragment& f)=0;
  /**
   * A member function that resets the parser and rewinds to the beginning of
   * the input.
   */
  virtual void reset() = 0;
};

/**
 * The Writer class is an abstract class for implementing a SAMWriter or
 * BAMWriter. It writes Fragment objects back to file (in SAM/BAM format) with
 * per-mapping probabilistic assignments, or by sampling a single mapping based
 * on assignment probabilities.
 **/
class Writer {
 protected:
  /**
   * A private bool that specifies if a single alignment should be sampled
   * (true) or all output with their respective posterior probabilities (false).
   */
  bool _sample;

 public:
  /**
   * Dummy destructor.
   */
  virtual ~Writer(){}
  /**
   * A member function that writes all mappings of the fragment to the ouptut
   * file along with their posterior probabilities in the "XP" field.
   * @param f the processed Fragment to output.
   */
  virtual void write_fragment(Fragment& f)=0;
};

/**
 * The BAMParser class fills Fragment objects by parsing an input file in BAM
 * format.
 **/
class BAMParser : public Parser {
  /**
   * A private pointer to the BamReader object which directly parses the BAM
   * file. Automatically deleted with BAMParser object.
   */
  boost::scoped_ptr<BamTools::BamReader> _reader;
  /**
   * A private member function to parse a single read alignment and store the
   * data in _read_buff.
   * @param alignment a BamAlignment containing the data parsed by BamTools.
   * @return True if the mapping is valid and false otherwise
   */
  bool map_end_from_alignment(BamTools::BamAlignment& alignment);

 public:
  /**
   * BAMParser constructor sets the reader.
   * @param reader a pointer to the BamReader object that will directly parse
   *        the BAM file.
   */
  BAMParser(BamTools::BamReader* reader);
  /**
   * An accessor for the header string.
   * @return The header string.
   */
  const std::string header() const { return _reader->GetHeaderText(); }
  /**
   * A member function that loads all mappings of the next fragment into the
   * given Fragment object.
   * @param f the empty Fragment to add mappings to.
   * @return True iff more reads remain in the BAM file/stream.
   */
  bool next_fragment(Fragment& f);
  /**
   * A member function that resets the parser and rewinds to the beginning of
   * the BAM file.
   */
  void reset();
};

/**
 * The SAMParser class fills Fragment objects by parsing an input in SAM format.
 * The input must come from a file.
 **/
class SAMParser : public Parser
{
  /**
   * A private pointer to the input stream in SAM format.
   */
  std::istream* _in;
  /**
   * A private string storing the SAM header.
   */
  std::string _header;
  /**
   * A private member function to parse a single read alignment and store the
   * data in _read_buff.
   * @param alignment a BamAlignment containing the data parsed by BamTools.
   * @return True if the mapping is valid and false otherwise
   */
  bool map_end_from_line(char* line);

public:
  /**
   * SAMParser constructor removes the header and parses the first line to
   * start the first Fragment.
   * @param in the input stream in SAM format, which must be a file
   */
  SAMParser(std::istream* in);
  /**
   * An accessor for the header string.
   * @return The header string.
   */
  const std::string header() const { return _header; }
  /**
   * A member function that loads all mappings of the next fragment into the
   * given Fragment object.
   * @param f the empty Fragment to add mappings to.
   * @return True iff more reads remain in the SAM/BAM file/stream.
   */
  bool next_fragment(Fragment& f);
  /**
   * A member function that resets the parser and rewinds to the beginning of
   * the SAM file.
   */
  void reset();
};

/**
 * The BAMWriter class writes Fragment objects back to file in BAM format with
 * per-mapping probabilistic assignments, or by sampling a single mapping based
 * on assignment probabilities.
 **/
class BAMWriter : public Writer {
  /**
   * A private pointer to the BamTools::BamWriter object which directly writes
   * the BAM file. Automatically deleted with BAMWriter object.
   */
  boost::scoped_ptr<BamTools::BamWriter> _writer;

 public:
  /**
   * BAMWriter constructor stores a pointer to the BamTools::BamWriter object
   * that will directly write to the BAM file.
   * @param writer pointer to the BamTools::BamWriter objected assocaited with
   *        the output BAM file.
   * @param sample specifies if a single alignment should be sampled based on
   *        posteriors (true) or all output with their respective posterior

   *        probabilities (false).
   */
  BAMWriter(BamTools::BamWriter* writer, bool sample);
  /**
   * BAMWriter destructor closes the BamTools::BamWriter object.
   */
  ~BAMWriter();
  /**
   * A member function that writes the mappings to the output BAM file. If
   * _sample is true, a only one alignment is output, otherwise all mappings are
   * output along with their probabilities in the "XP" field.
   * @param f the processed Fragment to output alignments of.
   */
  void write_fragment(Fragment& f);
};

/**
 * The SAMWriter class writes Fragment objects back to file in SAM format with
 * per-mapping probabilistic assignments, or by sampling a single mapping based
 * on assignment probabilities.

 **/
class SAMWriter : public Writer {
  /**
   * A private pointer to the output stream to which the alignments are
   * written in SAM format. Deleted with the SAMWriter object.
   */
  boost::scoped_ptr<std::ostream> _out;

 public:
  /**
   * SAMWriter constructor stores a pointer to the output stream.
   * @param out SAM output stream
   * @param sample specifies if a single alignment should be sampled based on
   *        posteriors (true) or all output with their respective posterior
   *        probabilities (false).
   */
  SAMWriter(std::ostream* out, bool sample);
  /**
   * SAMWriter destructor flushes the output stream.
   */
  ~SAMWriter();
  /**
   * A member function that writes the mappings to the output SAM file. If
   * _sample is true, a only one alignment is output, otherwise all mappings are
   * output along with their probabilities in the "XP" field.
   * @param f the processed Fragment to output alignments of.
   */
  void write_fragment(Fragment& f);
};

/**
 * The MapParser class is meant to be run as a separate thread from the main
 * processing. Once started, this thread will read input from a file or stream
 * in SAM/BAM format, parse, and collect read alignments into fragment
 * alignments, and fragment alignments into fragments, which are placed on a
 * buffer for the processing thread. Once the processing thread copies the
 * fragment address from the buffer, the parser is unlocked to load the next
 * fragment.  The process stops when EOF is reached.
 **/
class MapParser
{
  /**
   * A private pointer to the Parser object that will read the input in SAM/BAM
   * format. Automatically deleted with MapParser.
   */
  boost::scoped_ptr<Parser> _parser;
  /**
   * A private pointer to the Writer object that will write the output in
   * SAM/BAM format. Automatically deleted with MapParser.
   */
  boost::scoped_ptr<Writer> _writer;
  /**
   * A private pointer to other variables associated with the input.
   */
  Library* _lib;
  /**
   * A private boolean specifying whether to output the modified Fragments after
   * processing.
   */
  bool _write_active;

 public:
  /**
   * MapParser constructor determines what format the input is in and
   * initializes the correct parser and writer (if appropriate).
   * @param lib pointer to variables associated with the input, including file
   *        path.
   * @param write_active bool to initialize _write_active.
   */
  MapParser(Library* lib, bool write_active);
  /**
   * A member function that drives the parse thread. When all valid mappings of
   * a fragment have been parsed, its mapped targets are found and the
   * information is passed in a Fragment object to the processing thread through
   * a queue in the ParseThreadSafety struct. After processing, the Fragment
   * returns on a different in queue, and is written to the output map file
   * (depending on settings) and deleted.
   * @param thread_safety a pointer to the struct containing shared queues with
   *        the processing thread.
   * @param stop_at a size_t indicating how many reads to process before
   *        stopping (disabled if 0, default).
   */
  void threaded_parse(ParseThreadSafety* thread_safety, size_t stop_at=0);
  /**
   * A member function that counts the number of distinct fragments in the input file.
   * @return Number of fragments in the input file.
   */
  size_t count_fragments(size_t stop_at);
  /**
   * An accessor for the target name to index map. Returns a reference that does
   * not outlive this.
   * @return Reference to the target-to-index map.
   */
  const TargIndex& targ_index() { return _parser->targ_index(); }
  /**
   * An accessor for the set of target IDs with hits. Returns a reference that does
   * not outlive this.
   * @return Reference to the hit target ID set.
   */
  TargIDSet& hit_target_ids() { return _parser->hit_target_ids(); }
  /**
   * A member function that adds a new TargID to the set of targets with hits from
   * parsed fragmentsl.
   */
  void add_hit_target(TargID id) { (void) _parser->hit_target_ids().insert(id); }
  /**
   * An accessor for the target-to-length map. Returns a reference that does not
   * outlive this.
   * @return Reference to the target-to-length map.
   */
  const TargIndex& targ_lengths() { return _parser->targ_lengths(); }
  /**
   * A mutator for the write-active status of the parser. This specifies whether
   * or not the alignments (sampled or with probs) should be ouptut.
   * @param b updated write-active status
   */
  void write_active(bool b) { _write_active = b; }
  /**
   * A member function that resets the input parser.
   */
  void reset_reader() { _parser->reset(); }
};

#endif
