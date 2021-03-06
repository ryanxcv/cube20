#include "phase2prune.h"
#include <cstdio>
#include <iostream>
using namespace std;

struct corner_reduce {
  unsigned char m, parity;
  lookup_type c, minbits;
} corner_reduction[FACT8];
lookup_type edgeud_remap[KOCSYMM][FACT8];

int phase2prune::cornermax;
unsigned int phase2prune::memsize;
unsigned int *phase2prune::mem;

#ifdef HALF
const char *const phase2prune::filename = "p2p1h.dat";
#endif
#ifdef QUARTER
const char *const phase2prune::filename = "p2p1q.dat";
#endif
#ifdef SLICE
const char *const phase2prune::filename = "p2p1s.dat";
#endif
int phase2prune::file_checksum;

inline int corner_coordinate(const permcube &pc) {
  return (pc.c8_4 * FACT4 + pc.ctp) * FACT4 + pc.cbp;
}
inline int edge_coordinate(const permcube &pc) {
  return (permcube::c12_8[pc.et] * FACT4 + pc.etp) * FACT4 + pc.ebp;
}

static int datahash(unsigned int *dat, int sz, int seed) {
  while (sz > 0) {
    sz -= 4;
    seed = 37 * seed + *dat++;
  }
  return seed;
}

int phase2prune::lookup(const cubepos &cp) {
  permcube pc(cp);
  return lookup(pc);
}
int phase2prune::lookup(const permcube &pc) {
  int cc = corner_coordinate(pc);
  corner_reduce &cr = corner_reduction[cc];
  int off = cr.c * FACT8 + edgeud_remap[cr.m][edge_coordinate(pc)];
  int r = (mem[off >> 3] >> (4 * (off & 7))) & 0xf;
#ifdef QUARTER
  return cr.parity + 2 * r;
#else
  if (r == 0 && pc == identity_pc)
    return 0;
  else
    return r + 1;
#endif
}
int phase2prune::getindex(const permcube &pc) {
  int cc = corner_coordinate(pc);
  corner_reduce &cr = corner_reduction[cc];
  int off = cr.c * FACT8 + edgeud_remap[cr.m][edge_coordinate(pc)];
  return 2 * off + cr.parity;
}

void phase2prune::gen_table() {
  memset(mem, 255, memsize);
  cout << "Gen phase2" << flush;
#ifdef QUARTER
  mem[0] &= ~15;
  int seen = 1;
  for (int d = 1; d < 31; d++) {
    int backwards = (d >= 27);
    unsigned int seek = (d - 1) >> 1;
#else
  mem[0] &= ~14;
  int seen = 1;
  duration();
  for (int d = 0; d < 15; d++) {
    int backwards = (d >= 13);
    unsigned int seek = (d ? d - 1 : 1);
    int newval = d;
#endif
    for (int c8_4 = 0; c8_4 < C8_4; c8_4++)
      for (int ctp = 0; ctp < FACT4; ctp++)
        for (int cbp = 0; cbp < FACT4; cbp++) {
          permcube pc;
          pc.c8_4 = c8_4;
          pc.ctp = ctp;
          pc.cbp = cbp;
          int oc = corner_coordinate(pc);
          corner_reduce &cr = corner_reduction[oc];
#ifdef QUARTER
          if ((cr.minbits & 1) && (cr.parity != (d & 1))) {
#else
          if (cr.minbits & 1) {
#endif

            permcube pc2, pc3, pc4;
            cubepos cp2, cp3;
            int off = corner_reduction[oc].c * (FACT8 / 8);
            for (int mv = 0; mv < NMOVES_EXT; mv++) {
              if (!kocsymm::in_Kociemba_group(mv))
                continue;
#ifdef QUARTER
              int newval = d;
              if (mv >= NMOVES)
                newval++;
              newval >>= 1;
#endif
              pc2 = pc;
              pc2.move(mv);
              int dest_off = corner_coordinate(pc2);
              corner_reduce &cr = corner_reduction[dest_off];
              int destat = cr.c * (FACT8 / 8);
              for (int m = cr.m; (1 << m) <= cr.minbits; m++)
                if ((cr.minbits >> m) & 1) {

                  int at = 0;
                  for (int e8_4 = 0; e8_4 < C8_4; e8_4++) {
                    int et = permcube::c8_12[e8_4];
                    int t1 = permcube::eperm_move[et][mv];
                    int eb =
                        kocsymm::epsymm_compress[0xf0f -
                                                 kocsymm::epsymm_expand[et]];
                    int t2 = permcube::eperm_move[eb][mv] & 31;
                    int dst1 = permcube::c12_8[t1 >> 5] * 24 * 24;
                    t1 &= 31;
                    if (backwards) {
                      for (int etp = 0; etp < FACT4; etp++)
                        for (int ebpo = 0; ebpo < FACT4; ebpo += 8) {
                          unsigned int v = mem[off + (at >> 3)];
                          v &= v >> 1;
                          if ((0x11111111 & v & (v >> 2)) != 0) {
                            for (int ebpi = 0; ebpi < 8; ebpi++, at++)
                              if (((mem[off + (at >> 3)] >> (4 * (at & 7))) &
                                   0xf) == 0xf) {
                                int ebp = ebpo + ebpi;

                                int etp1 = permcube::s4mul[etp][t1];
                                int ebp1 = permcube::s4mul[ebp][t2];
                                int dat =
                                    edgeud_remap[m][dst1 + etp1 * 24 + ebp1];
                                unsigned int val = (mem[destat + (dat >> 3)] >>
                                                    (4 * (dat & 7))) & 0xf;
                                if (val == seek) {
                                  mem[off + (at >> 3)] -= (0xf - newval)
                                                          << (4 * (at & 7));
                                  seen++;
                                }
                              }
                          } else {
                            at += 8;
                          }
                        }
                    } else {
                      for (int etp = 0; etp < FACT4; etp++)
                        for (int ebpo = 0; ebpo < FACT4; ebpo += 8) {
                          if (mem[off + (at >> 3)] != 0xffffffff) {
                            for (int ebpi = 0; ebpi < 8; ebpi++, at++)
                              if (((mem[off + (at >> 3)] >> (4 * (at & 7))) &
                                   0xf) == seek) {
                                int ebp = ebpo + ebpi;

                                int etp1 = permcube::s4mul[etp][t1];
                                int ebp1 = permcube::s4mul[ebp][t2];
                                int dat =
                                    edgeud_remap[m][dst1 + etp1 * 24 + ebp1];
                                int val = (mem[destat + (dat >> 3)] >>
                                           (4 * (dat & 7))) &
                                          0xf;
                                if (val == 0xf) {
                                  mem[destat + (dat >> 3)] -=
                                      (0xf - newval) << (4 * (dat & 7));
                                  seen++;
                                }
                              }
                          } else {
                            at += 8;
                          }
                        }
                    }
                  }
                }
            }
          }
        }
#ifndef QUARTER
    if (d == 0)
      mem[0] &= ~15;
#endif
    cout << " " << d << flush;
  }
  cout << " done." << endl << flush;
}

void phase2prune::write_table() {
  FILE *f = fopen(filename, "wb");
  if (f == 0)
    error("! cannot write pruning file to current directory");
  if (fwrite(mem, 1, memsize, f) != memsize)
    error("! error writing pruning table");
  if (fwrite(&file_checksum, sizeof(int), 1, f) != 1)
    error("! error writing pruning table");
  fclose(f);
}

void phase2prune::check_integrity() {
  if (file_checksum != datahash(mem, memsize, 0))
    error("! integrity of pruning table compromised");
  cout << "Verified integrity of phase two pruning data: " << file_checksum
       << endl;
}

moveseq phase2prune::solve(const permcube &pc, int maxlen) {
  moveseq r;
  for (int d = lookup(pc); d <= maxlen; d++)
    if (solve(pc, d, CANONSEQSTART, r)) {
      reverse(r.begin(), r.end());
      break;
    }
  return r;
}

int phase2prune::solve(const permcube &pc, int togo, int canonstate,
                       moveseq &r) {
  if (lookup(pc) > togo)
    return 0;
  if (pc == identity_pc)
    return 1;
  if (togo-- <= 0)
    return 0;
  permcube pc2;
#ifdef QUARTER
  int mask = cubepos::cs_mask_ext(canonstate) & 0xf0c3;
#else
  int mask = cubepos::cs_mask(canonstate) & 0227227227;
#endif
  while (mask) {
    int ntogo = togo;
    int mv = ffs1(mask);
    mask &= mask - 1;
#ifdef QUARTER
    if (mv >= NMOVES) {
      if (togo == 0)
        continue;
      ntogo = togo - 1;
    }
#endif
    pc2 = pc;
    pc2.move(mv);
    if (solve(pc2, ntogo, cubepos::next_cs(canonstate, mv), r)) {
#ifdef QUARTER
      if (mv >= NMOVES) {
        int nmv = mv - NMOVES;
        nmv = 2 * (nmv + 1 + nmv / 2);
        r.push_back(nmv);
        r.push_back(nmv);
      } else {
        r.push_back(mv);
      }
#else
      r.push_back(mv);
#endif
      return 1;
    }
  }
  return 0;
}

#include "table.h"

void phase2prune::init(int suppress_writing) {
  static int initialized = 0;
  if (initialized)
    return;
  initialized = 1;

  cornermax = 0;
  for (int c8_4 = 0; c8_4 < C8_4; c8_4++)
    for (int ctp = 0; ctp < FACT4; ctp++)
      for (int cbp = 0; cbp < FACT4; cbp++) {
        permcube pc;
        pc.c8_4 = c8_4;
        pc.ctp = ctp;
        pc.cbp = cbp;
        int oc = corner_coordinate(pc);
        int minc = oc;
        int minm = 0;
        int minbits = 1;
        cubepos cp;
        pc.set_perm(cp);
        for (int m = 1; m < 16; m++) {
          cubepos cp2;
          cp.remap_into(m, cp2);
          permcube pc2(cp2);
          int tc = corner_coordinate(pc2);
          if (tc < minc) {
            minc = tc;
            minm = m;
            minbits = 1 << m;
          } else if (tc == minc)
            minbits |= 1 << m;
        }
        corner_reduce &cr = corner_reduction[oc];
        if (oc == minc)
          cr.c = cornermax++;
        cr.m = minm;
        cr.c = corner_reduction[minc].c;
        cr.minbits = minbits;
        cr.parity = (permcube::c8_4_parity[c8_4] + ctp + cbp) & 1;
      }

  int at = 0;
  cubepos cp, cp2;
  for (int e8_4 = 0; e8_4 < C8_4; e8_4++) {
    permcube pc;
    pc.et = permcube::c8_12[e8_4];
    pc.eb = kocsymm::epsymm_compress[0xf0f - kocsymm::epsymm_expand[pc.et]];
    for (int etp = 0; etp < FACT4; etp++) {
      pc.etp = etp;
      for (int ebp = 0; ebp < FACT4; ebp++, at++) {
        pc.ebp = ebp;
        for (int m = 0; m < KOCSYMM; m++) {
          pc.set_edge_perm(cp);
          cp.remap_into(m, cp2);
          permcube pc2(cp2);
          edgeud_remap[m][at] = edge_coordinate(pc2);
        }
      }
    }
  }

  memsize = cornermax * (FACT8 / 2);
  read_table_req((unsigned char **) &mem, filename, memsize);
}
