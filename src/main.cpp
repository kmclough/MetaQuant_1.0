/**
 *  main.cpp
 *  MetaQuant
 *
 *  Created by Kevin McLoughlin on 8/16/14
 *  Based on eXpress, created by Adam Roberts on 3/23/11.
 *  Copyright 2014 Kevin McLoughlin, Adam Roberts. All rights reserved.
 **/

#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <boost/thread.hpp>
#include <boost/chrono/chrono.hpp>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "main.h"
#include "bundles.h"
#include "targets.h"
#include "lengthdistribution.h"
#include "fragments.h"
#include "biascorrection.h"
#include "mismatchmodel.h"
#include "mapparser.h"
#include "threadsafety.h"
#include "robertsfilter.h"
#include "directiondetector.h"
#include "library.h"

#ifndef WIN32
  #include "update_check.h"
#endif

using namespace std;
namespace po = boost::program_options;
namespace fs = boost::filesystem;

Logger logger;

// the forgetting factor parameter controls the growth of the fragment mass
double ff_param = 0.85;
// The learning rate parameters serve a similar function for the online variational
// Bayes updates
double learn_rate_offset = 1000.0;
double learn_rate_power = 0.85;

// Alpha and beta parameters for the stick-breaking Beta prior
double alpha_prior = 1.0;
double beta_prior = 1.0;

// the burn-in parameter determines how many reads are required before the
// error and bias models are applied to probabilistic assignment
size_t burn_in = 50000;
size_t burn_out = 5000000;
bool burned_out = false;

// =-= This will undoubtedly have to increase in the future
size_t max_read_len = 250;

size_t max_indel_size = 10;

size_t stop_at = 0;
size_t trace_frags = 0;

// file location parameters
string output_dir = ".";
string fasta_file_name = "";
string in_map_file_names = "";
string param_file_name = "";

// intial pseudo-count parameters (non-logged)
double expr_alpha = .005;
double fld_alpha = 1;
double bias_alpha = 1;
double mm_alpha = 1;

size_t bias_model_order = 3;

// fragment length parameters
size_t def_fl_max = 800;
size_t def_fl_mean = 200;
size_t def_fl_stddev = 80;
size_t def_fl_kernel_n = 4;
double def_fl_kernel_p = 0.5;

// option parameters
bool edit_detect = false;
bool error_model = true;
bool bias_correct = true;
bool calc_covar = false;
bool output_align_prob = false;
bool output_align_samp = false;
bool output_running_rounds = false;
bool output_running_reads = false;
size_t num_threads = 2;
size_t library_size = 0;
int bias_thread_sleep_time = 300;

// directional parameters
Direction direction = BOTH;

bool running = true;

// used for multiple rounds of VBayes
bool first_round = true;
bool last_round = true;
bool batch_mode = false;
bool online_additional = false;
bool both = false;
size_t remaining_rounds = 0;

bool spark_pre = false;


/**
 * Parses argument options and sets variables appropriately.
 * @param ac number of arguments.
 * @param pointer to array of arguments as character arrays.
 * @return True iff there was an error.
 */
bool parse_options(int ac, char ** av) {

  size_t additional_online = 0;
  size_t additional_batch = 0;
  
  po::options_description standard("Standard Options");
  standard.add_options()
  ("help,h", "produce help message")
  ("output-dir,o", po::value<string>(&output_dir)->default_value(output_dir),
   "write all output files to this directory")
  ("frag-len-mean,m", po::value<size_t>(&def_fl_mean)->default_value(def_fl_mean),
   "prior estimate for average fragment length")
  ("frag-len-stddev,s",
   po::value<size_t>(&def_fl_stddev)->default_value(def_fl_stddev),
   "prior estimate for fragment length std deviation")
  ("additional-batch,B",
   po::value<size_t>(&additional_batch)->default_value(additional_batch),
   "number of additional batch VBayes rounds after initial online round")
  ("additional-online,O",
   po::value<size_t>(&additional_online)->default_value(additional_online),
   "number of additional online VBayes rounds after initial online round")
  ("max-read-len,L",
   po::value<size_t>(&max_read_len)->default_value(max_read_len),
   "maximum allowed length of a read")
  ("output-align-prob",
   "output alignments (sam/bam) with probabilistic assignments")
  ("output-align-samp",
   "output alignments (sam/bam) with sampled assignments")
  ("fr-stranded",
   "accept only forward->reverse alignments (second-stranded protocols)")
  ("rf-stranded",
   "accept only reverse->forward alignments (first-stranded protocols)")
  ("f-stranded",
   "accept only forward single-end alignments (second-stranded protocols)")
  ("r-stranded",
   "accept only reverse single-end alignments (first-stranded protocols)")
  // =-= TODO: If we get a dedicated support web site for MetaQuant, uncomment the
  // =-= following line.
  //("no-update-check", "disables automatic check for update via web")
  ("logtostderr", "prints all logging messages to stderr")
  ;
  
  po::options_description advanced("Advanced Options");
  advanced.add_options()
  ("forget-param,f", po::value<double>(&ff_param)->default_value(ff_param),
   "sets the 'forgetting factor' parameter (0.5 < c <= 1)")
  ("learn-power-param", po::value<double>(&learn_rate_power)->default_value(learn_rate_power),
   "sets the learning rate exponent parameter (0.5 < power <= 1)")
  ("learn-offset-param", po::value<double>(&learn_rate_offset)->default_value(learn_rate_offset),
   "sets the learning rate offset parameter (> 0)")
  ("library-size", po::value<size_t>(&library_size),
   "specifies library size for FPKM instead of calculating from alignments")
  ("max-indel-size",
   po::value<size_t>(&max_indel_size)->default_value(max_indel_size),
   "sets the maximum allowed indel size, affecting geometric indel prior")
  ("calc-covar", "calculate and output covariance matrix")
  ("expr-alpha", po::value<double>(&expr_alpha)->default_value(expr_alpha),
   "sets the strength of the prior, per bp")
  ("trace-frags", po::value<size_t>(&trace_frags)->default_value(trace_frags),
   "sets the interval between trace messages in terms of fragments processed, disabled with 0")
  ("stop-at", po::value<size_t>(&stop_at)->default_value(stop_at),
   "sets the number of fragments to process, disabled with 0")
  ("alpha-prior", po::value<double>(&alpha_prior)->default_value(alpha_prior),
   "sets alpha parameter for prior in stick-breaking Beta distribution")
  ("beta-prior", po::value<double>(&beta_prior)->default_value(beta_prior),
   "sets beta parameter for prior in stick-breaking Beta distribution")
  ("burn-in", po::value<size_t>(&burn_in)->default_value(burn_in),
   "sets number of fragments after which to begin using auxiliary parameters")
  ("burn-out", po::value<size_t>(&burn_out)->default_value(burn_out),
   "sets number of fragments after which to stop updating auxiliary parameters")
  ("bias-thread-sleep-time", po::value<int>(&bias_thread_sleep_time)->default_value(bias_thread_sleep_time),
   "sets time interval in seconds to wait between executions of the bias update thread")
  ("no-bias-correct", "disables bias correction")
  ("no-error-model", "disables error modelling")
  ("aux-param-file",
   po::value<string>(&param_file_name)->default_value(param_file_name),
   "path to file containing auxiliary parameters to use instead of learning")
  ;


  po::options_description hidden("Experimental/Debug Options");
  hidden.add_options()
  ("num-threads,p", po::value<size_t>(&num_threads)->default_value(num_threads),
   "number of threads (>= 2)")
  ("edit-detect","")
  ("single-round", "")
  ("output-running-rounds", "")
  ("output-running-reads", "")
  ("batch-mode","")
  ("both","")
  ("sam-file", po::value<string>(&in_map_file_names)->default_value(""), "")
  ("fasta-file", po::value<string>(&fasta_file_name)->default_value(""), "")
  ("bias-model-order",
   po::value<size_t>(&bias_model_order)->default_value(bias_model_order),
   "sets the order of the Markov chain used to model sequence bias")
  ;

  po::positional_options_description positional;
  positional.add("fasta-file",1).add("sam-file",1);

  po::options_description cmdline_options;
  cmdline_options.add(standard).add(advanced).add(hidden);

  bool error = false;
  po::variables_map vm;
  try {
    po::store(po::command_line_parser(ac, av).options(cmdline_options)
              .positional(positional).run(), vm);
  } catch (po::error& e) {
    logger.info("Command-Line Argument Error: %s.", e.what());
    error = true;
  }
  po::notify(vm);

  if (ff_param > 1.0 || ff_param < 0.5) {
    logger.info("Command-Line Argument Error: forget-param/f option must be "
                "between 0.5 and 1.0.");
    error= true;
  }

  if (learn_rate_power > 1.0 || learn_rate_power < 0.5) {
    logger.info("Command-Line Argument Error: learn_power-param option must be "
                "between 0.5 and 1.0.");
    error= true;
  }

  if (learn_rate_offset <= 0.0) {
    logger.info("Command-Line Argument Error: learn-offset-param option must be "
                "greater than 0.");
    error= true;
  }

  if (fasta_file_name == "") {
    logger.info("Command-Line Argument Error: target sequence fasta file "
                "required.");
    error = true;
  }

  if (error || vm.count("help")) {
    cerr << "metaquant v" << PACKAGE_VERSION << endl
         << "-----------------------------\n"
         << "Usage:  metaquant [options] <target_seqs.fa> <hits.(sam/bam)>\n"
         << "Required arguments:\n"
         << " <target_seqs.fa>     target sequence file in fasta format\n"
         << " <hits.(sam/bam)>     read alignment file in SAM or BAM format\n\n"
         << standard
         << advanced;
    return 1;
  }

  if (param_file_name.size()) {
    burn_in = 0;
    burn_out = 0;
    burned_out = true;
  }
  
  size_t stranded_count = 0;
  if (vm.count("fr-stranded")) {
    direction = FR;
    stranded_count++;
  }
  if (vm.count("rf-stranded")) {
    direction = RF;
    stranded_count++;
  }
  if (vm.count("f-stranded")) {
    direction = F;
    stranded_count++;
  }
  if (vm.count("r-stranded")) {
    direction = R;
    stranded_count++;
  }
  if (stranded_count > 1) {
    logger.severe("Multiple strandedness flags cannot be specified in the same "
                  "run.");
  }
  if (vm.count("logtostderr")) {
    logger.info_out(&cerr);
  }
  
  edit_detect = vm.count("edit-detect");
  calc_covar = vm.count("calc-covar");
  bias_correct = !(vm.count("no-bias-correct"));
  error_model = !(vm.count("no-error-model"));
  output_align_prob = vm.count("output-align-prob");
  output_align_samp = vm.count("output-align-samp");
  output_running_rounds = vm.count("output-running-rounds");
  output_running_reads = vm.count("output-running-reads");
  batch_mode = vm.count("batch-mode");
  both = vm.count("both");
  remaining_rounds = max(additional_online, additional_batch);
  spark_pre = vm.count("preprocess");

  if (batch_mode) {
    // =-= TODO Implement batch mode
    logger.severe("Batch mode is not yet supported.");
    ff_param = 1;
  }
  if (additional_batch > 0) {
    // =-= TODO Implement batch mode
    logger.severe("Batch mode is not yet supported.");
  }
  
  if (additional_online > 0 && additional_batch > 0) {
    logger.severe("Cannot add both online and batch rounds.");
  } else if (additional_online > 0) {
    online_additional = true;
  }
  
  if (output_align_prob && output_align_samp) {
    logger.severe("Cannot output both alignment probabilties and sampled "
                  "alignments.");
  }
  if ((output_align_prob || output_align_samp) && remaining_rounds == 0) {
    logger.warn("It is recommended that at least one additional round "
                "be used when outputting alignment probabilities or sampled "
                "alignments. Use the '-B' or '-O' option to enable.");
  }
  
  // We have 1 processing thread and 1 parsing thread always, so we should not
  // count these as additional threads.
  if (num_threads < 2) {
    num_threads = 0;
  }
  num_threads -= 2;
  if (num_threads > 0) {
    num_threads -= edit_detect;
  }
  if (in_map_file_names == "") {
    logger.severe("MetaQuant cannot process streaming input; you must provide a SAM file name.");
  }
  if (remaining_rounds) {
    last_round = false;
  }

#ifndef WIN32
  // =-= TODO: Uncomment this code and replace URLs in update_check.h, if we ever
  // =-= get a dedicated support web site for MetaQuant.
  //if (!vm.count("no-update-check")) {
  //  check_version(PACKAGE_VERSION);
  //}
#endif

  return 0;
}

/**
 * This function writes the current abundance parameters to one file and the
 * auxiliary parameters for each library to a separate file.
 * @param libs a Librarian containing the parameters tables for each Library.
 * @param tot_counts a size_t for the total number of fragments processed thus
          far.
 * @param n an int suffix to add to the output subdirectory. No subdirectory is
 *        used if -1 (default).
 */
void output_results(Librarian& libs, size_t tot_counts, int n=-1) {
  char buff[500];
  string dir = output_dir;
  if (n >= 0) {
    sprintf(buff, "%s/x_%d", output_dir.c_str(), n);
    logger.info("Writing results to %s.", buff);
    dir = string(buff);
    try {
      fs::create_directories(dir);
    } catch (fs::filesystem_error& e) {
      logger.severe(e.what());
    }
  }
  libs[0].targ_table->output_results(dir, tot_counts);

  for (size_t l = 0; l < libs.size(); l++) {
    if (libs.size() > 1) {
      sprintf(buff, "%s/params_metaquant_%d.txt", dir.c_str(), (int)l+1);
    } else {
      sprintf(buff, "%s/params_metaquant.txt", dir.c_str());
    }
    ofstream paramfile(buff);
    (libs[l].fld)->append_output(paramfile, "Fragment");
    if (libs[l].mismatch_table) {
      (libs[l].mismatch_table)->append_output(paramfile);
    }
    if (libs[l].bias_table) {
      (libs[l].bias_table)->append_output(paramfile);
    }
    paramfile.close();
  }
}

/**
 * This function handles the probabilistic assignment of multi-mapped reads. The
 * marginal likelihoods are calculated for each mapping, and the mass of the
 * fragment is divided based on the normalized marginals to update the model
 * parameters.
 * @param frag_p pointer to the fragment to probabilistically assign.
 */
void process_fragment(Fragment* frag_p) {
  Fragment& frag = *frag_p;
  const Library& lib = *frag.lib();

  //boost::chrono::high_resolution_clock::time_point pfrag_start = 
  //   boost::chrono::high_resolution_clock::now();

  // sort hits to avoid deadlock
  // =-= TODO: Is this necessary?
  frag.sort_hits();
  double mass_n = frag.mass();

  assert(frag.num_hits());

  size_t num_solvable = 0;

  // Update bundles and merge in first loop
  Bundle* bundle = frag.hits()[0]->target()->bundle();
  
  if (frag.num_hits() > 1) {
    // Calculate likelihoods for each alignment
    for (size_t i = 0; i < frag.num_hits(); ++i) {
      FragHit& hit = *frag.hits()[i];
      Target* t = hit.target();
      FragTarget* frag_targ = frag.get_frag_target(hit.target_id());
      
      bundle = lib.targ_table->merge_bundles(bundle, t->bundle());
      t->bundle(bundle);
      
      hit.params()->align_likelihood = hit.align_likelihood(frag, t);
      frag_targ->total_align_likelihood = log_add(frag_targ->total_align_likelihood,
                                                  hit.params()->align_likelihood);
      num_solvable += t->solvable();
    }
  } else {
    FragHit& hit = *frag.hits()[0];
    FragTarget* frag_targ = frag.get_frag_target(hit.target_id());
    frag_targ->total_align_likelihood = LOG_1;
  }

  if (first_round) {
    bundle->incr_counts();
  }
  if (first_round || online_additional) {
    bundle->incr_mass(mass_n);
  }
  
  // Perform variational Bayes updates for posterior probabilities that current fragment
  // came from each target, and stick-breaking beta distribution parameters for each target
  // hit by fragment
  //boost::chrono::high_resolution_clock::time_point vbayes_start = 
  //   boost::chrono::high_resolution_clock::now();
  frag.var_bayes_updates();
  //boost::chrono::high_resolution_clock::time_point vbayes_end = 
  //   boost::chrono::high_resolution_clock::now();
  //long vbayes_time = 
  //   (long)((vbayes_end-vbayes_start).count());

  // Update the posterior likelihoods for each specific alignment; use these to update
  // the parameters for the bias, mismatch and length distributions.
  double total_likelihood = LOG_0;

  for (size_t i = 0; i < frag.num_hits(); ++i) {
    FragHit& hit = *frag[i];
    FragTarget* frag_targ = frag.get_frag_target(hit.target_id());

    hit.params()->full_likelihood = hit.params()->align_likelihood 
                                     + frag_targ->expect_log_sample_prob;
    
    total_likelihood = log_add(total_likelihood, hit.params()->full_likelihood);
  }

  if (islzero(total_likelihood)){
    logger.warn("Fragment '%s' has 0 likelihood of originating from any target in the "
                "target database. Skipping...", frag.name().c_str());
    return;
  }

  // Normalize the fragment hit posteriors to add to 1
  for (size_t i = 0; i < frag.num_hits(); ++i) {
    FragHit& hit = *frag[i];
    Target* t  = hit.target();
    double p = hit.params()->full_likelihood-total_likelihood;
    hit.params()->posterior = p;

    // Update bias, mismatch, fld parameters
    if (first_round) {
      double r = rand()/double(RAND_MAX);
      
      if (i == 0 || frag[i-1]->target_id() != t->id()) {
        t->incr_counts(frag.num_frag_targets() <= 1);
      }
      if (!t->solvable() && num_solvable == frag.num_hits()-1) {
        t->solvable(true);
      }
      if (edit_detect && lib.mismatch_table) {
        (lib.mismatch_table)->update(hit, p, lib.mass_n);
      }
      if (!burned_out && r < sexp(p)) {
        if (lib.mismatch_table && !edit_detect) {
          (lib.mismatch_table)->update(hit, LOG_1, lib.mass_n);
        }
        if (hit.pair_status() == PAIRED) {
          (lib.fld)->add_val(hit.length(), lib.mass_n);
        }
        if (lib.bias_table) {
          (lib.bias_table)->update_observed(hit, lib.mass_n);
        }
      }
    }
    if (calc_covar && (last_round || online_additional)) {
      double var = 2*mass_n + p + log_sub(LOG_1, p);
      lib.targ_table->update_covar(hit.target_id(), hit.target_id(), var);
      for (size_t j = i+1; j < frag.num_hits(); ++j) {
        const FragHit& m2 = *frag.hits()[j];
        double p2 = m2.params()->full_likelihood-total_likelihood;
        if (sexp(p2) == 0) {
          continue;
        }
        double covar = 2*mass_n + p + p2;
        lib.targ_table->update_covar(hit.target_id(), m2.target_id(), covar);
      }
    }
  }
  //boost::chrono::high_resolution_clock::time_point pfrag_end = 
  //   boost::chrono::high_resolution_clock::now();
  //long pfrag_time = 
  //   (long)((pfrag_end-pfrag_start).count()) - vbayes_time;
  //logger.info("Fragment %d: %d hits; var_bayes_update() time = %ld ns", 
  //   lib.n, frag.num_hits(), vbayes_time);
  //logger.info("  remaining process_fragment() time = %ld ns", pfrag_time);
}

/**
 * This function processes Fragments asynchronously. Fragments are popped from
 * a threadsafe input queue, processed, and then pushed onto a threadsafe output
 * queue.
 * @param pts pointer to a struct with the input and output Fragment queues.
 */
void proc_thread(ParseThreadSafety* pts) {
  while (true) {
    Fragment* frag = pts->proc_on.pop();
    if (!frag) {
      break;
    }
    process_fragment(frag);
    pts->proc_out.push(frag);
  }
}

/**
 * This is the driver function for the main processing thread. This function
 * updates the current fragment mass for libraries, dispatches fragments to be
 * processed once they are passed by the parsing thread, outputs intermediate
 * results, and handles additional online rounds.
 * @param libs a struct containing pointers to the parameter tables (bias_table,
 *        mismatch_table, fld) and parser for all libraries being processed.
 * @return The total number of fragments processed.
 */
size_t threaded_calc_abundances(Librarian& libs) {
  logger.info("Processing input fragment alignments...");
  boost::scoped_ptr<boost::thread> bias_update;

  size_t n = 1;
  size_t num_frags = 0;
  double mass_n = 0;
  double learning_rate = pow(learn_rate_offset, -learn_rate_power);

  // For log-scale output
  size_t i = 1;
  size_t j = 6;

  DirectionDetector dir_detector;
  Fragment* frag;
  
  while (true) {
    // Loop through libraries
    for (size_t l = 0; l < libs.size(); l++) {
      Library& lib = libs[l];
      libs.set_curr(l);
      MapParser& map_parser = *lib.map_parser;
      boost::mutex bu_mut;
      // Used to signal bias update thread
      running = true;
      ParseThreadSafety pts(max((int)num_threads,10));
      boost::thread parse(&MapParser::threaded_parse, &map_parser, &pts, stop_at);
      vector<boost::thread*> thread_pool;
      RobertsFilter frags_seen;

      burned_out = lib.n >= burn_out;
      while(true) {
        if (lib.n == burn_in) {
          bias_update.reset(new boost::thread(&TargetTable::asynch_bias_update,
                                              lib.targ_table, &bu_mut, 
                                              bias_thread_sleep_time));
          if (lib.mismatch_table) {
            (lib.mismatch_table)->activate();
          }
        }
        if (lib.n == burn_out) {
          if (lib.mismatch_table) {
            (lib.mismatch_table)->fix();
          };
          burned_out = true;
        }
        // Start threads once aux parameters are burned out
        if (burned_out && num_threads && thread_pool.size() == 0) {
          lib.targ_table->enable_bundle_threadsafety();
          thread_pool = vector<boost::thread*>(num_threads);
          for (size_t k = 0; k < thread_pool.size(); k++) {
            thread_pool[k] = new boost::thread(proc_thread, &pts);
          }
        }

        // Pop next parsed fragment and set its forgetting mass and learning rate
        frag = pts.proc_in.pop();
        if (frag) {
          frag->mass(mass_n);
          frag->learning_rate(learning_rate);
          dir_detector.add_fragment(frag);
        }

        // Test that we have not already seen this fragment
        if (frag && first_round && frags_seen.test_and_push(frag->name())) {
          logger.severe("Alignments are not properly sorted. Read '%s' has "
                        "alignments which are non-consecutive.",
                        frag->name().c_str());
        }

        // If multi-threaded and burned out, push to the processing queue
        if (num_threads && burned_out) {
          // If no more fragments, send stop signal (NULL) to processing threads
          if (!frag) {
            for (size_t k = 0; k < thread_pool.size(); ++k) {
              pts.proc_on.push(NULL);
            }
            break;
          }
          pts.proc_on.push(frag);
        } else {
          if (!frag) {
            break;
          }
          {
            // Block the bias update thread from updating the parameter tables
            // during processing. We don't need to do this during multi-threaded
            // processing since the parameters are burned out before we start
            // the threads.
            boost::unique_lock<boost::mutex> lock(bu_mut);
            process_fragment(frag);
            pts.proc_out.push(frag);
          }
        }

        // Output intermediate results, if necessary
        if (output_running_reads && n == i*pow(10.,(double)j)) {
          boost::unique_lock<boost::mutex> lock(bu_mut);
          lib.targ_table->estimate_abundances();
          output_results(libs, n, (int)n);
          if (i++ == 9) {
            i = 1;
            j++;
          }
        }
        num_frags++;

        // Output progress
        if ((trace_frags > 0) && (num_frags % trace_frags == 0)) {
          logger.info("Fragments Processed (%s): %d\tNumber of Bundles: %d.",
                      lib.in_file_name.c_str(), num_frags,
                      lib.targ_table->num_bundles());
          dir_detector.report_if_improper_direction();
        }

        n++;
        lib.n++;
        mass_n += ff_param*log((double)n-1) - log(pow(n,ff_param) - 1);
        lib.mass_n += ff_param*log((double)lib.n-1) -
                      log(pow(lib.n,ff_param) - 1);
        learning_rate = pow(learn_rate_offset + n, -learn_rate_power);
      }

      // Signal bias update thread to stop
      running = false;

      parse.join();
      foreach(boost::thread* t, thread_pool) {
        t->join();
      }

      lib.targ_table->disable_bundle_threadsafety();
      lib.targ_table->collapse_bundles();
      
      if (bias_update) {
        logger.info("Waiting for auxiliary parameter update to complete...");
        bias_update->join();
        bias_update.reset(NULL);
      }
    }

    if (online_additional && remaining_rounds--) {
      if (output_running_rounds) {
        libs[0].targ_table->estimate_abundances();
        output_results(libs, n, (int)remaining_rounds);
      }

      logger.info("%d remaining rounds.", remaining_rounds);
      first_round = false;
      last_round = (remaining_rounds==0 && !both);
      for (size_t l = 0; l < libs.size(); l++) {
        libs[l].map_parser->write_active(last_round);
        libs[l].map_parser->reset_reader();
      }
      num_frags = 0;
    } else {
      break;
    }
  }

  logger.info("COMPLETED: Processed %d mapped fragments, targets are in %d "
              "bundles.", num_frags, libs[0].targ_table->num_bundles());

  return num_frags;
}

/**
 * The main function instantiates the library parameter tables and parsers,
 * calls the processing function, and outputs the results. Also handles
 * additional batch rounds.
 */
int estimation_main() {
  
  if (output_dir != ".") {
    try {
      fs::create_directories(output_dir);
    } catch (fs::filesystem_error& e) {
      logger.info(e.what());
    }
  }
  
  if (!fs::exists(output_dir)) {
    logger.severe("Cannot create directory %s.", output_dir.c_str());
  }
  
  // Parse input file names and instantiate Library structs.
  vector<string> file_names;
  char buff[999];
  strcpy(buff, in_map_file_names.c_str());
  char * pch = strtok (buff,",");
  while (pch != NULL) {
    file_names.push_back(pch);
    pch = strtok (NULL, ",");
  }
  if (file_names.size() == 0) {
    file_names.push_back("");
  }
  Librarian libs(file_names.size());
  for (size_t i = 0; i < file_names.size(); ++i) {
    char out_map_file_name[500] = "";
    if (output_align_prob) {
      sprintf(out_map_file_name, "%s/hits.%d.prob",
              output_dir.c_str(), (int)i+1);
    }
    if (output_align_samp) {
      sprintf(out_map_file_name, "%s/hits.%d.samp",
              output_dir.c_str(), (int)i+1);
    }
    
    libs[i].in_file_name = file_names[i];
    libs[i].out_file_name = out_map_file_name;
    libs[i].map_parser.reset(new MapParser(&libs[i], last_round));
    libs[i].total_frags = libs[i].map_parser->count_fragments(stop_at);
    libs[i].map_parser->reset_reader();

    if (param_file_name.size()) {
      libs[i].fld.reset(new LengthDistribution(param_file_name, "Fragment"));
      libs[i].mismatch_table.reset((error_model) ?
                                              new MismatchTable(param_file_name)
                                              : NULL);
      libs[i].bias_table.reset((bias_correct) ? new BiasBoss(bias_model_order,
                                                         param_file_name):NULL);
    } else {
      libs[i].fld.reset(new LengthDistribution(fld_alpha, def_fl_max,
                                               def_fl_mean, def_fl_stddev,
                                               def_fl_kernel_n,
                                               def_fl_kernel_p));
      libs[i].mismatch_table.reset((error_model) ? new MismatchTable(mm_alpha)
                                                   :NULL);
      libs[i].bias_table.reset((bias_correct) ? new BiasBoss(bias_model_order,
                                                    bias_alpha):NULL);
    }
    if (i > 0 &&
        (libs[i].map_parser->targ_index() != libs[i-1].map_parser->targ_index()
         || libs[i].map_parser->targ_lengths() !=
         libs[i-1].map_parser->targ_lengths())) {
          logger.severe("Alignment file headers do not match for '%s' and '%s'.",
                        file_names[i-1].c_str(), file_names[i].c_str());
        }
  }
  
  boost::shared_ptr<TargetTable> targ_table(
                                  new TargetTable(fasta_file_name,
                                                  edit_detect,
                                                  param_file_name.size(),
                                                  &libs));
  size_t max_target_length = 0;
  for(size_t tid=0; tid < targ_table->size(); tid++) {
    max_target_length = max(max_target_length,
                            targ_table->get_targ(tid)->length());
  }

  for (size_t i = 0; i < libs.size(); ++i) {
    libs[i].targ_table = targ_table;
    if (bias_correct) {
      libs[i].bias_table->copy_expectations(*(libs.curr_lib().bias_table));
    }
  }
  double num_targ = (double)targ_table->size();
  
  if (calc_covar && (double)SSIZE_MAX < num_targ*(num_targ+1)) {
    logger.warn("Your system is unable to represent large enough values for "
                "efficiently hashing target pairs. Covariance calculation will "
                "be disabled.");
    calc_covar = false;
  }
  
  if (batch_mode) {
    targ_table->round_reset();
  }
  
  // Execute the main abundance calculation loop

  size_t tot_counts = threaded_calc_abundances(libs);
  if (library_size) {
    tot_counts = library_size;
  }
  
  if (!burned_out && bias_correct && param_file_name == "") {
    logger.warn("Not enough fragments observed to accurately learn bias "
                "parameters. Either disable bias correction "
                "(--no-bias-correct) or provide a file containing auxiliary "
                "parameters (--aux-param-file).");
  }
  
  if (both) {
    remaining_rounds = 1;
    online_additional = false;
  }
  
  targ_table->estimate_abundances();
  // =-= TODO: Add support for additional batch VBayes rounds
  //targ_table->round_reset();
  ff_param = 1.0;

  first_round = false;
  
  // =-= TODO: The following loop never gets executed, I think, since remaining_rounds
  // =-= should already be 0 on exiting threaded_calc_abundances(). Get rid of it?
  while (!last_round) {
    remaining_rounds--;
    logger.info("\nRe-estimating counts with additional round of VBayes (%d "
                "remaining)...", remaining_rounds);
    last_round = (remaining_rounds == 0);
    for (size_t l = 0; l < libs.size(); l++) {
      libs[l].map_parser->write_active(last_round);
      libs[l].map_parser->reset_reader();
    }
    tot_counts = threaded_calc_abundances(libs);
    if (library_size) {
      tot_counts = library_size;
    }
    targ_table->estimate_abundances();
    if (output_running_rounds) {
      output_results(libs, tot_counts, (int)remaining_rounds);
    }
    targ_table->round_reset();
  }
  
  logger.info("Writing results to file...");
  targ_table->estimate_abundances();
  output_results(libs, tot_counts);
  logger.info("Done.");
  
  return 0;
}


int main (int argc, char ** argv)
{

  srand((unsigned int)time(NULL));
  int parse_ret = parse_options(argc,argv);
  if (parse_ret) {
    return parse_ret;
  }
  
  return estimation_main();
}
