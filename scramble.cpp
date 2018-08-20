#include "twophase.h"

// Generate a random integer from 0 (inclusive) to n (exclusive)
// TODO: better RNG, e.g. urandom
int randint(int n) {
    return rand() % n;
}

/*
 * Translated from Shuang Chen's min2phase
 * https://github.com/cs0x7f/min2phase
 */

int getnparity(int idx, int n) {
  int p = 0;
  for (int i = 2; i <= n + 1; i++) {
    p ^= idx % i;
    idx /= i;
  }
  return p & 1;
}

void setnperm(uint8_t* arr, int idx, int n) {
  long val = 0xFEDCBA9876543210L;
  long extract = 0;
  for (int p = 2; p <= n; p++) {
    extract = extract << 4 | idx % p;
    idx /= p;
  }
  for (int i = 0; i < n - 1; i++) {
    int v = ((int) extract & 0xf) << 2;
    extract >>= 4;
    arr[i] = val >> v & 0xf;
    long m = (1L << v) - 1;
    val = (val & m) | (val >> 4 & ~m);
  }
  arr[n - 1] = val & 0xf;
}

void setcperm(uint8_t* ca, int idx) {
    setnperm(ca, idx, 8);
}

void settwist(uint8_t* co, int twist) {
  int i, parity = 0;
  for (i = 6; i >= 0; i--) {
      parity += co[i] = (char) (twist % 3);
      twist /= 3;
  }
  co[7] = (char) ((3 - parity % 3) % 3);
}

void setflip(uint8_t* eo, short flip) {
    int i, parity = 0;
    for (i = 10; i >= 0; i--) {
        parity += eo[i] = (char) (flip % 2);
        flip /= 2;
    }
    eo[11] = (char) ((2 - parity % 2) % 2);
}

cubepos randomstate() {
  time_t t;
  srand((unsigned) time(&t));

  int epval, cpval = randint(40320);
  int parity = getnparity(cpval, 8);
  do {
      epval = randint(479001600);
  } while (getnparity(epval, 12) != parity);

  cubepos c;
  unsigned char cp[8],  co[8];
  unsigned char ep[12], eo[12];

  setcperm(cp, cpval);
  setnperm(ep, epval, 12);
  settwist(co, randint(2187));
  setflip( eo, randint(2048));

  for (int i = 0; i < 8; i++)
    c.c[i] = cubepos::corner_val(cp[i], co[i]);
  for (int i = 0; i < 12; i++)
    c.e[i] = cubepos::edge_val(ep[i], eo[i]);
  return c;
}

void pretty(int mv) {
  char f;
  char *p = &f;
  cubepos::append_face(p, mv / TWISTS);
  cout << f;
  int t = mv % TWISTS;
  if (t)
#ifndef QUARTER
      cout << " 2'"[t];
#else
      cout << " '"[t];
#endif
}

void prettyprint_moveseq(const moveseq &seq) {
    pretty(seq[0]);
    for (unsigned int i = 1; i < seq.size(); i++) {
        cout << " ";
        pretty(seq[i]);
    }
    cout << endl;
}

int main(int argc, char *argv[]) {
  double progstart = walltime();
  duration();

  twophasesolver solver;

  const int skipwrite = 0;
  phase1prune::init(skipwrite);
  phase2prune::init(skipwrite);

  cubepos cp = randomstate();
  static int input_seq = 1;
  int seq = input_seq++;
  solution s = solver.solve(seq, cp);
  display(s);
  prettyprint_moveseq(s.moves);

  cerr << "Completed in " << (walltime() - progstart) << endl;
}
