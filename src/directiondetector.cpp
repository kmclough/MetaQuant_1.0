/*
 *  directiondetector.cpp
 *  MetaQuant
 *
 *  Kevin McLoughlin
 *  Based on code from eXpress, created by Adam Roberts in 2012.
 *  Copyright 2014 Kevin McLoughlin, Adam Roberts. All rights reserved.
 */

#include "directiondetector.h"

#include <iostream>
#include "fragments.h"
#include "main.h"

using namespace std;

DirectionDetector::DirectionDetector()
: _num_fr(0), _num_rf(0), _num_f(0), _num_r(0) {}

void DirectionDetector::add_fragment(Fragment* f) {
  foreach (FragHit* h, f->hits()) {
    switch(h->pair_status()) {
      case PAIRED: {
        if (h->left_read() == h->first_read()) {
          _num_fr++;
        } else {
          assert(h->right_read() == h->first_read());
          _num_rf++;
        }
        break;
      }
      case LEFT_ONLY: {
        _num_f++;
        break;
      }
      case RIGHT_ONLY: {
        _num_r++;
        break;
      }
    }
  }
}

bool DirectionDetector::report_if_improper_direction() {
  size_t num_paired = _num_fr + _num_rf;
  size_t num_single = _num_f + _num_r;
  if (num_paired + num_single == 0) {
    return false;
  }
  if (num_paired == 0) {
    // Single-end case
    double max_dir = max(_num_f, _num_r);
    double min_dir = min(_num_f, _num_r);
    if (min_dir < max_dir / 2) {
      if (_num_f > _num_r && direction != F) {
        logger.warn("The observed alignments appear disporportionately on "
                    "the forward strand (%d  vs. %d). If your library is "
                    "strand-specific and single-end, you should use the "
                    "--f-stranded option to avoid incorrect results.",
                    _num_f, _num_r);
        return true;
      } else if (_num_f < _num_r && direction != R) {
        logger.warn("The observed alignments appear disporportionately on "
                    "the reverse strand (%d vs. %d). If your library is "
                    "strand-specific and single-end, you should use the "
                    "--r-stranded option to avoid incorrect results.",
                    _num_r, _num_f);
        return true;
      }
    }
  } else {
    // Paired-end case
    size_t fr = _num_f + _num_fr;
    size_t rf = _num_r + _num_rf;
    double max_dir = max(fr, rf);
    double min_dir = min(fr, rf);
    if (min_dir < max_dir / 2) {
      if (fr > rf && direction != FR) {
        logger.warn("The observed alignments appear disporportionately in the"
                    "the forward-reverse order (%d vs %d). If your library is "
                    "strand-specific, you should use the --fr-stranded option "
                    "to avoid incorrect results.", fr, rf);
        return true;
      } else if (rf > fr && direction != RF) {
        logger.warn("The observed alignments appear disporportionately in "
                    "the reverse-forward order (%d vs. %d). If your library is "
                    "strand-specific, you should use the --rf-stranded option "
                    "to avoid incorrect results.", rf, fr);
        return true;
      }
    }
  }
  return false;
}


