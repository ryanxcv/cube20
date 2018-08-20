#pragma once

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <stddef.h>
#include <vector>
#ifdef _WIN32
#include <intrin.h>
#include <windows.h>
inline int ffs1(int v) {
  unsigned long r;
  _BitScanForward(&r, v);
  return (int)r;
}
#else
inline int ffs1(int v) { return ffs(v) - 1; }
#include <sys/time.h>
#endif
using namespace std;

#ifdef HALF
#ifdef SLICE
#error "Can't define HALF and SLICE"
#endif
#ifdef QUARTER
#error "Can't define HALF and SLICE"
#endif
#ifdef AXIAL
#error "Can't define HALF and AXIAL"
#endif
#else
#ifdef SLICE
#ifdef QUARTER
#error "Can't define SLICE and QUARTER"
#endif
#ifdef AXIAL
#error "Can't define SLICE and AXIAL"
#endif
#else
#ifdef QUARTER
#ifdef AXIAL
#error "Can't define SLICE and AXIAL"
#endif
#else
#ifndef AXIAL
#error "Please define one of HALF, SLICE, QUARTER, or AXIAL"
#endif
#endif
#endif
#endif

#ifdef HALF
const int NMOVES = 18;
const int TWISTS = 3;
#endif
#ifdef QUARTER
const int NMOVES = 12;
const int TWISTS = 2;
#endif
#ifdef SLICE
const int NMOVES = 27;
const int TWISTS = 3;
#endif
#ifdef AXIAL
const int NMOVES = 45;
const int TWISTS = 3;
#endif
const int FACES = 6;
const int M = 48;
const int CUBIES = 24;

extern const class cubepos identity_cube;

#ifdef QUARTER
const int NMOVES_EXT = NMOVES + 4;
#else
const int NMOVES_EXT = NMOVES;
#endif

typedef vector<int> moveseq;

const long long ALLMOVEMASK = (1LL << NMOVES) - 1;
const long long ALLMOVEMASK_EXT = (1LL << NMOVES_EXT) - 1;

#ifdef HALF
const int CANONSEQSTATES = FACES + 1;
#endif
#ifdef QUARTER
const int CANONSEQSTATES = 2 * FACES + 1;
#endif
#ifdef SLICE
const int CANONSEQSTATES = 5 * FACES / 2 + 1;
#endif
#ifdef AXIAL
const int CANONSEQSTATES = 3 + 1;
#endif
const int CANONSEQSTART = 0;

void error(const char *s);
double myrand();
inline int random_move() { return (int)(NMOVES * myrand()); }
inline int random_move_ext() { return (int)(NMOVES_EXT * myrand()); }
double walltime();
double duration();

void init_mutex();
void get_global_lock();
void release_global_lock();
#ifdef THREADS
#ifdef _WIN32
#include <process.h>
#include <windows.h>
#define THREAD_RETURN_TYPE unsigned int
#define THREAD_DECLARATOR __stdcall
#else
#include <pthread.h>
#define THREAD_RETURN_TYPE void *
#define THREAD_DECLARATOR
#endif
const int MAX_THREADS = 128;
void spawn_thread(int i, THREAD_RETURN_TYPE(THREAD_DECLARATOR *p)(void *),
                  void *o);
void join_thread(int i);
#else
#define THREAD_RETURN_TYPE void *
#define THREAD_DECLARATOR
const int MAX_THREADS = 1;
#endif

class cubepos {
public:
  inline bool operator<(const cubepos &cp) const {
    return memcmp(this, &cp, sizeof(cp)) < 0;
  }
  inline bool operator==(const cubepos &cp) const {
    return memcmp(this, &cp, sizeof(cp)) == 0;
  }
  inline bool operator!=(const cubepos &cp) const {
    return memcmp(this, &cp, sizeof(cp)) != 0;
  }

  static inline int edge_perm(int cubieval) { return cubieval >> 1; }
  static inline int edge_ori(int cubieval) { return cubieval & 1; }
  static inline int corner_perm(int cubieval) { return cubieval & 7; }
  static inline int corner_ori(int cubieval) { return cubieval >> 3; }
  static inline int edge_flip(int cubieval) { return cubieval ^ 1; }
  static inline int edge_val(int perm, int ori) { return perm * 2 + ori; }
  static inline int corner_val(int perm, int ori) { return ori * 8 + perm; }
  static inline int edge_ori_add(int cv1, int cv2) {
    return cv1 ^ edge_ori(cv2);
  }
  static inline int corner_ori_add(int cv1, int cv2) {
    return mod24[cv1 + (cv2 & 0x18)];
  }
  static inline int corner_ori_sub(int cv1, int cv2) {
    return cv1 + corner_ori_neg_strip[cv2];
  }
  static void init();

  inline cubepos(const cubepos &cp = identity_cube) { *this = cp; }
  cubepos(int, int, int);

  void move(int mov);

  static int invert_move(int mv) { return inv_move[mv]; }
  static moveseq invert_sequence(const moveseq &sequence);
  void invert_into(cubepos &dst) const;

  void movepc(int mov);

  static void mul(const cubepos &a, const cubepos &b, cubepos &r);
  inline static void mulpc(const cubepos &a, const cubepos &b, cubepos &r) {
    mul(b, a, r);
  }

  static void skip_whitespace(const char *&p);
  static int parse_face(const char *&p);
  static int parse_face(char f);
#if defined(SLICE) || defined(AXIAL)
  static int parse_moveface(const char *&p);
  static int parse_moveface(char f);
  static void append_moveface(char *&p, int f) { *p++ = movefaces[f]; }
  static void append_face(char *&p, int f) { *p++ = movefaces[f]; }
#else
  static void append_face(char *&p, int f) { *p++ = faces[f]; }
#endif
  static int parse_move(const char *&p);
  static void append_move(char *&p, int mv);
  static moveseq parse_moveseq(const char *&p);
  static void append_moveseq(char *&p, const moveseq &seq);
  static char *moveseq_string(const moveseq &seq);

  const char *parse_Singmaster(const char *p);
  char *Singmaster_string() const;

  void remap_into(int m, cubepos &dst) const;
  void canon_into48(cubepos &dst) const;
  void canon_into48_aux(cubepos &dst) const;
  void canon_into96(cubepos &dst) const;

  void randomize();

  static inline int next_cs(int cs, int mv) { return canon_seq[cs][mv]; }
  static inline long long cs_mask(int cs) { return canon_seq_mask[cs]; }
  static inline long long cs_mask_ext(int cs) { return canon_seq_mask_ext[cs]; }

  static unsigned char corner_ori_inc[CUBIES], corner_ori_dec[CUBIES],
      corner_ori_neg_strip[CUBIES], mod24[2 * CUBIES];

  static char faces[FACES];
#ifdef SLICE
  static char movefaces[FACES + 3];
#endif
#ifdef AXIAL
  static char movefaces[FACES + 9];
#endif

  static unsigned char edge_trans[NMOVES_EXT][CUBIES],
      corner_trans[NMOVES_EXT][CUBIES];

  static unsigned char inv_move[NMOVES_EXT];

  static unsigned char face_map[M][FACES], move_map[M][NMOVES_EXT];
  static unsigned char invm[M], mm[M][M];
  static unsigned char rot_edge[M][CUBIES], rot_corner[M][CUBIES];

  static unsigned char canon_seq[CANONSEQSTATES][NMOVES_EXT];
  static long long canon_seq_mask[CANONSEQSTATES];
  static long long canon_seq_mask_ext[CANONSEQSTATES];

  unsigned char c[8];

  unsigned char e[12];
};

static cubepos cubepos_initialization_hack(1, 2, 3);
