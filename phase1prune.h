#pragma once

#include "kocsymm.h"

#ifdef HALF
const int BYTES_PER_ENTRY = 4;
#endif
#ifdef SLICE
const int BYTES_PER_ENTRY = 5;
#endif
#ifdef QUARTER
const int BYTES_PER_ENTRY = 4;
#endif

class phase1prune {
public:
  static void init(int suppress_writing = 0);
  static int lookup(const kocsymm &kc, int &mask);

  static void gen_table();
  static int read_table();
  static void write_table();
  static void check_integrity();

  static int lookup(const kocsymm &kc);
  static int lookup(const kocsymm &kc, int togo, int &nextmovemask);
  static moveseq solve(kocsymm kc);

  static unsigned int memsize;
  static unsigned char *mem;
  static int file_checksum;
  static const char *const filename;
};
