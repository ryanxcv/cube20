#include "phase1prune.h"
#include <cstdio>
#include <iostream>
using namespace std;

unsigned int phase1prune::memsize;
unsigned char *phase1prune::mem;
int phase1prune::file_checksum;
#ifdef HALF
const char *const phase1prune::filename = "p1p1h.dat";
#endif
#ifdef QUARTER
const char *const phase1prune::filename = "p1p1q.dat";
#endif
#ifdef SLICE
const char *const phase1prune::filename = "p1p1s.dat";
#endif

#ifdef SLICE
static unsigned char map_phase1_offsets[KOCSYMM][6];
#else
static unsigned char map_phase1_offsets[KOCSYMM][3];
#endif
static int map_phase1[2][12][256];

static int datahash(unsigned int *dat, int sz, int seed) {
  while (sz > 0) {
    sz -= 4;
    seed = 37 * seed + *dat++;
  }
  return seed;
}

void phase1prune::gen_table() {
  memset(mem, -1, memsize);
  mem[0] = 0;
  int seen = 1;
  cout << "Gen phase1" << flush;
  for (int d = 1;; d++) {
    int lastiter = (seen == CORNERRSYMM * EDGEOSYMM * EDGEPERM);
    int seek = d - 1;
    int at = 0;
    for (int cs = 0; cs < CORNERRSYMM; cs++) {
      int csymm = kocsymm::cornersymm_expand[cs];
      for (int eosymm = 0; eosymm < EDGEOSYMM; eosymm++)
        for (int epsymm = 0; epsymm < EDGEPERM; epsymm++, at += BYTES_PER_ENTRY)
#ifdef SLICE
          if ((mem[at] >> 4) == seek) {
#else
          if (mem[at] == seek) {
#endif
            int deltadist[NMOVES];
            for (int mv = 0; mv < NMOVES; mv++) {
              int rd = 0;
              kocsymm kc(csymm, eosymm, epsymm);
              kc.move(mv);
              corner_mapinfo &cm = kocsymm::cornersymm[kc.csymm];
              for (int m = cm.minmap; cm.minbits >> m; m++)
                if ((cm.minbits >> m) & 1) {
                  int deosymm =
                      kocsymm::edgeomap[kocsymm::edgepxor[kc.epsymm][m >> 3] ^
                                        kc.eosymm][m];
                  int depsymm = kocsymm::edgepmap[kc.epsymm][m];
                  int dat =
                      ((cm.csymm * EDGEOSYMM + deosymm) * EDGEPERM + depsymm) *
                      BYTES_PER_ENTRY;
                  rd = mem[dat];
                  if (rd == 255) {
                    rd = d;
#ifdef SLICE
                    mem[dat] = rd << 4;
#else
                    mem[dat] = rd;
#endif
                    seen++;
#ifdef SLICE
                  } else {
                    rd >>= 4;
#endif
                  }
                }
              deltadist[mv] = rd - seek;
            }
#ifdef SLICE
            mem[at] = mem[at] & 0xf0;
            mem[at + 1] = 0;
#endif
            for (int b = 0; b < 3; b++) {
              int v = 0;
#ifdef SLICE
              int clim = 2;
#else
              int clim = 1;
#endif
              for (int c = clim; c >= 0; c--) {
                int vv = 0;
#ifdef QUARTER
                for (int t = 1; t >= 0; t--)
                  vv = 3 * vv + deltadist[2 * b + 6 * c + t] + 1;
                v = 9 * v + vv;
#else
                int cnts[3];
                cnts[0] = cnts[1] = cnts[2] = 0;
                for (int t = 2; t >= 0; t--) {
                  vv = 2 * vv + deltadist[3 * b + 9 * c + t];
                  cnts[1 + deltadist[3 * b + 9 * c + t]]++;
                }
                if (cnts[0] > 0 && cnts[2] > 0) {
                  cout << "counts are " << cnts[0] << " " << cnts[1] << " "
                       << cnts[2] << endl;
                  error("! bad delta distance values within one face turn set");
                }
                if (cnts[0])
                  vv += 7;
                else
                  vv += 8;
                v = 16 * v + vv;
#endif
              }
#ifdef SLICE
              mem[at + b + 2] = v;
              mem[at + (b + 1) / 2] |= (v >> 8) << (4 * (b & 1));
#else
              mem[at + b + 1] = v;
#endif
            }
          }
    }
    cout << " " << d << flush;
    if (lastiter)
      break;
  }
  cout << " done." << endl << flush;
}

void phase1prune::write_table() {
  FILE *f = fopen(filename, "wb");
  if (f == 0)
    error("! cannot write pruning file to current directory");
  if (fwrite(mem, 1, memsize, f) != memsize)
    error("! error writing pruning table");
  if (fwrite(&file_checksum, sizeof(int), 1, f) != 1)
    error("! error writing pruning table");
  fclose(f);
}

void phase1prune::check_integrity() {
  if (file_checksum != datahash((unsigned int *)mem, memsize, 0))
    error("! integrity of pruning table compromised");
  cout << "Verified integrity of phase one pruning data: " << file_checksum
       << endl;
}

int phase1prune::lookup(const kocsymm &kc) {
  corner_mapinfo &cm = kocsymm::cornersymm[kc.csymm];
  int m = cm.minmap;
  int r = mem[BYTES_PER_ENTRY *
              (((cm.csymm * EDGEOSYMM) +
                kocsymm::edgeomap[kocsymm::edgepxor[kc.epsymm][m >> 3] ^
                                  kc.eosymm][m]) *
                   495 +
               kocsymm::edgepmap[kc.epsymm][m])];
#ifdef SLICE
  return r >> 4;
#else
  return r;
#endif
}

int phase1prune::lookup(const kocsymm &kc, int togo, int &nextmovemask) {
  corner_mapinfo &cm = kocsymm::cornersymm[kc.csymm];
  int m = cm.minmap;
  int off = BYTES_PER_ENTRY *
            (((cm.csymm * EDGEOSYMM) +
              kocsymm::edgeomap[kocsymm::edgepxor[kc.epsymm][m >> 3] ^
                                kc.eosymm][m]) *
                 495 +
             kocsymm::edgepmap[kc.epsymm][m]);
  int r = mem[off];
#ifdef SLICE
  r >>= 4;
#endif
  if (togo < r) {
    nextmovemask = 0;
  } else if (togo > r + 1) {
    nextmovemask = ALLMOVEMASK;
  } else {
    int(*p)[256] = map_phase1[togo - r];
    unsigned char *o = map_phase1_offsets[m];
#ifdef SLICE
    nextmovemask = p[o[0]][mem[off + 2]] + p[o[1]][mem[off + 3]] +
                   p[o[2]][mem[off + 4]] +
                   (((p[o[3]][mem[off] & 15] + p[o[4]][mem[off + 1] >> 4] +
                      p[o[5]][mem[off + 1] & 15]) &
                     0777)
                    << 18);
#else
    nextmovemask =
        p[o[0]][mem[off + 1]] + p[o[1]][mem[off + 2]] + p[o[2]][mem[off + 3]];
#endif
  }
  return r;
}

#include "table.h"

void phase1prune::init(int suppress_writing) {
  static int initialized = 0;
  if (initialized)
    return;
  initialized = 1;

  memsize = BYTES_PER_ENTRY * CORNERRSYMM * EDGEOSYMM * EDGEPERM;
  read_table_req(&mem, filename, memsize);

  for (int m = 0; m < KOCSYMM; m++) {
    for (int f = 0; f < 3; f++) {
      int mv = f * TWISTS;
      int mv2 = cubepos::move_map[m][mv];
      int f2 = mv2 / TWISTS;
      int key = 0;
      if (mv2 % TWISTS == TWISTS - 1)
        key++;
      if (f2 >= 3)
        key += 2;
      key += 4 * (f2 % 3);
      map_phase1_offsets[cubepos::invm[m]][f] = key;
#ifdef SLICE
      map_phase1_offsets[cubepos::invm[m]][f + 3] =
          (key & 0xd) ^ ((key & 2) >> 1);
#endif
    }
  }

  for (int slack = 0; slack < 2; slack++) {
    for (int key = 0; key < 12; key++) {
#ifdef QUARTER
      int nv[9];
      for (int nyb = 0; nyb < 9; nyb++) {
        int bits = 0;
        if (nyb % 3 <= slack)
          bits |= 1;
        if (nyb / 3 <= slack)
          bits |= 2;
        if (key & 1)
          bits = ((bits + 4 * bits) >> 1) & 3;
        if (key & 2)
          bits <<= 3 * TWISTS;
        bits <<= TWISTS * (key >> 2);
        nv[nyb] = bits;
      }
      int *a = map_phase1[slack][key];
      for (int byte = 0; byte < 81; byte++)
        a[byte] =
            nv[byte % 9] |
            (((nv[byte / 9] << (3 * TWISTS)) | (nv[byte / 9] >> (3 * TWISTS))) &
             ALLMOVEMASK);
#else
      int nv[16];
      for (int nyb = 0; nyb < 16; nyb++) {
        int bits = 0;
        if (slack && nyb <= 7) {
          bits = 7;
        } else if (slack == 0 && nyb >= 7) {
          bits = 0;
        } else {
          bits = 7 - (nyb & 7);
        }
        if (key & 1)
          bits = ((bits & 1) << 2) + (bits & 2) + (bits >> 2);
        if (key & 2)
          bits <<= 3 * TWISTS;
        bits <<= TWISTS * (key >> 2);
        nv[nyb] = bits;
      }
      int *a = map_phase1[slack][key];
      for (int byte = 0; byte < 256; byte++)
        a[byte] = nv[byte & 15] | (((nv[byte >> 4] << (3 * TWISTS)) |
                                    (nv[byte >> 4] >> (3 * TWISTS))) &
                                   0777777);
#endif
    }
  }
}

moveseq phase1prune::solve(kocsymm kc) {
  moveseq r;
  int d = phase1prune::lookup(kc);
  while (d > 0) {
    int nmm = 0;
    int t = phase1prune::lookup(kc, d, nmm);
    if (t == 0)
      break;
    if (t != d)
      error("! did not make progress");
    if (nmm == 0)
      error("! no solution?");
    int mv = ffs1(nmm);
    r.push_back(mv);
    kc.move(mv);
    d--;
  }
  return r;
}
