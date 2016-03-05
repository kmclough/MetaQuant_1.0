/*
 * memdebug.h
 *
 * Platform dependent function to determine amount of memory allocated, for debugging
 * memory problems
 *
 *  MetaQuant
 *
 *  Kevin McLoughlin
 *  Copyright 2015 Kevin McLoughlin. All rights reserved.
 */

#ifndef Linux
// MacOS (Darwin) version
#include <malloc/malloc.h>

size_t get_mem_used() {
  struct mstats stats = mstats();
  return stats.bytes_used;
}

size_t get_mem_free() {
  struct mstats stats = mstats();
  return stats.bytes_free;
}

#else // Linux

#include <malloc.h>

size_t get_mem_used() {
  struct mallinfo minfo = mallinfo();
  return (size_t) minfo.uordblks;
}

size_t get_mem_free() {
  struct mallinfo minfo = mallinfo();
  return (size_t) minfo.fordblks;
}

#endif // Linux
