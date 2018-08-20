#ifndef KOCSYMM_H
#define KOCSYMM_H
#include "cubepos.h"

const int CORNERSYMM = 2187;
const int EDGEOSYMM = 2048;
const int EDGEPERM = 495;
const int KOCSYMM = 16;
const int CORNERRSYMM = 168;

struct corner_mapinfo {
  unsigned short minbits;
  unsigned char csymm, minmap;
};

typedef unsigned short lookup_type;
class kocsymm {
public:
  kocsymm() : csymm(0), eosymm(0), epsymm(0) {}
  kocsymm(int c, int eo, int ep) : csymm(c), eosymm(eo), epsymm(ep) {}
  kocsymm(int) : csymm(0), eosymm(0), epsymm(0) { init(); }
  static void init();

  inline bool operator<(const kocsymm &kc) const {
    if (csymm != kc.csymm)
      return csymm < kc.csymm;
    if (eosymm != kc.eosymm)
      return eosymm < kc.eosymm;
    return epsymm < kc.epsymm;
  }
  inline bool operator==(const kocsymm &kc) const {
    return kc.csymm == csymm && kc.eosymm == eosymm && kc.epsymm == epsymm;
  }
  inline bool operator!=(const kocsymm &kc) const {
    return kc.csymm != csymm || kc.eosymm != eosymm || kc.epsymm != epsymm;
  }

  void move(int mv) {
    csymm = cornermove[csymm][mv];
    eosymm = edgeomove[eosymm][mv];
    epsymm = edgepmove[epsymm][mv];
  }

  kocsymm(const cubepos &cp);
  void set_coset(cubepos &cp);

  void canon_into(kocsymm &kc) const;

  int calc_symm() const;

  static inline int in_Kociemba_group(int mv) { return edgepmove[0][mv] == 0; }

  static lookup_type cornermove[CORNERSYMM][NMOVES_EXT];
  static lookup_type edgeomove[EDGEOSYMM][NMOVES_EXT];
  static lookup_type edgepmove[EDGEPERM][NMOVES_EXT];

  static lookup_type epsymm_compress[1 << 12];
  static lookup_type epsymm_expand[EDGEOSYMM];

  static lookup_type cornersymm_expand[CORNERRSYMM];
  static corner_mapinfo cornersymm[CORNERSYMM];
  static lookup_type edgeomap[EDGEOSYMM][KOCSYMM];
  static lookup_type edgepmap[EDGEPERM][KOCSYMM];
  static lookup_type edgepxor[EDGEPERM][2];

  lookup_type csymm, eosymm, epsymm;
};

static kocsymm identity_kc(1);

const int FACT4 = 24;
const int C8_4 = 70;
class permcube {
public:
  permcube();

  inline bool operator<(const permcube &pc) const {
    if (et != pc.et)
      return et < pc.et;
    if (em != pc.em)
      return em < pc.em;
    if (eb != pc.eb)
      return eb < pc.eb;
    if (etp != pc.etp)
      return etp < pc.etp;
    if (emp != pc.emp)
      return emp < pc.emp;
    if (ebp != pc.ebp)
      return ebp < pc.ebp;
    if (c8_4 != pc.c8_4)
      return c8_4 < pc.c8_4;
    if (ctp != pc.ctp)
      return ctp < pc.ctp;
    return cbp < pc.cbp;
  }
  inline bool operator==(const permcube &pc) const {
    return et == pc.et && em == pc.em && eb == pc.eb && etp == pc.etp &&
           emp == pc.emp && ebp == pc.ebp && c8_4 == pc.c8_4 && ctp == pc.ctp &&
           cbp == pc.cbp;
  }
  inline bool operator!=(const permcube &pc) const {
    return et != pc.et || em != pc.em || eb != pc.eb || etp != pc.etp ||
           emp != pc.emp || ebp != pc.ebp || c8_4 != pc.c8_4 || ctp != pc.ctp ||
           cbp != pc.cbp;
  }

  void move(int mv);

  void init_edge_from_cp(const cubepos &cp);
  void init_corner_from_cp(const cubepos &cp);
  permcube(const cubepos &cp);
  void set_edge_perm(cubepos &cp) const;
  void set_corner_perm(cubepos &cp) const;
  void set_perm(cubepos &cp) const;

  static void init();
  static unsigned char s4inv[FACT4];
  static unsigned char s4mul[FACT4][FACT4];
  static unsigned char s4compress[256];
  static unsigned char s4expand[FACT4];

  static unsigned char c8_4_compact[256];
  static unsigned char c8_4_expand[C8_4];
  static unsigned char c8_4_parity[C8_4];

  static unsigned char c12_8[EDGEPERM];
  static lookup_type c8_12[C8_4];

  static unsigned short eperm_move[EDGEPERM][NMOVES_EXT];
  static int cperm_move[C8_4][NMOVES_EXT];

  unsigned short et, em, eb;
  unsigned char etp, emp, ebp;
  unsigned char c8_4, ctp, cbp;
};

static permcube identity_pc;

#endif
