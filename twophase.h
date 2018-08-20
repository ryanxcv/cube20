#include "phase1prune.h"
#include "phase2prune.h"
#include <cstdio>
#include <iostream>
#include <map>
using namespace std;

int verbose = 1;
int numthreads = 1;
int numsols = 1;

int target_length = 20;
long long phase2limit = 0xffffffffffffffLL;
long long phase2total = 0LL;

int skipwrite = 0;

int axesmask = 63;

class solution {
public:
  solution(const cubepos &cparg, int seqarg, long long p2parg,
           moveseq &solarg) {
    cp = cparg;
    seq = seqarg;
    phase2probes = p2parg;
    moves = solarg;
  }
  solution() {}
  cubepos cp;
  int seq;
  long long phase2probes;
  moveseq moves;
};
map<int, solution> queue;
int next_sequence = 1;
int missed_target = 0;
int solved = 0;

#ifdef QUARTER
const int MAX_MOVES = 46;
#else
const int MAX_MOVES = 32;
#endif

int sloweq(const cubepos &cp1, const cubepos &cp2) {
  cubepos cp3;
  for (int m = 0; m < KOCSYMM; m++) {
    cp2.remap_into(m, cp3);
    if (cp1 == cp3)
      return 1;
  }
  return 0;
}

void display(const cubepos &cp, int seq, long long phase2probes, moveseq sol) {
  phase2total += phase2probes;
  if (verbose || (int)sol.size() > target_length) {
    if ((int)sol.size() > target_length)
      cout << "WARNING: missed target for " << cp.Singmaster_string() << endl;
    cout << "Solution " << seq << " len " << sol.size() << " probes "
         << phase2probes << endl;
    cout << cubepos::moveseq_string(sol) << endl;
  }
}

void display(solution& s) {
  display(s.cp, s.seq, s.phase2probes, s.moves);
}

class twophasesolver {
public:
  twophasesolver() {}
  cubepos pos;
  long long phase2probes;
  int bestsol;
  int keepbound;
  int keepcounts[MAX_MOVES];
  int keepsum;
  int finished;
  int curm;
  int solmap;
  int seq;

  unsigned char moves[MAX_MOVES];
  unsigned char bestmoves[MAX_MOVES];

  kocsymm kc6[6], kccanon6[6];
  cubepos cp6[6];
  permcube pc6[6];
  int mindepth[6];
  char uniq[6];
  int minmindepth;

  solution solve(int seqarg, cubepos &cp) {
    pos = cp;
    phase2probes = 0;
    bestsol = MAX_MOVES;
    keepbound = MAX_MOVES - 1;
    finished = 0;
    seq = seqarg;
    if (numsols > 1) {
      for (int i = 0; i < MAX_MOVES; i++)
        keepcounts[i] = 0;
      keepsum = 0;
    }

    minmindepth = MAX_MOVES;
    cubepos cpi, cp2;
    pos.invert_into(cpi);
    int ind = 0;
    for (int inv = 0; inv < 2; inv++)
      for (int mm = 0; mm < 3; mm++, ind++) {
        int m = KOCSYMM * mm;
        if (inv)
          cpi.remap_into(m, cp2);
        else
          pos.remap_into(m, cp2);
        cp6[ind] = cp2;
        kc6[ind] = kocsymm(cp2);
        pc6[ind] = permcube(cp2);
        kc6[ind].canon_into(kccanon6[ind]);
        mindepth[ind] = phase1prune::lookup(kc6[ind]);
        if (mindepth[ind] < minmindepth)
          minmindepth = mindepth[ind];
        uniq[ind] = 1 & (axesmask >> ind);
        for (int i = 0; i < ind; i++)
          if (uniq[i] && kccanon6[ind] == kccanon6[i] &&
              sloweq(cp6[ind], cp6[i])) {
            uniq[ind] = 0;
            break;
          }
        if (verbose > 1) {
          get_global_lock();
          cout << "Axis " << ind << " depth " << mindepth[ind] << " uniq "
               << (int)uniq[ind] << endl;
          release_global_lock();
        }
      }

    for (int d = minmindepth; d <= keepbound && !finished; d++) {
      for (curm = 0; curm < 6; curm++)
        if (uniq[curm] && d <= keepbound && !finished && d >= mindepth[curm]) {
          if (verbose > 1) {
            get_global_lock();
            cout << "Orientation " << curm << " at depth " << d << endl;
            release_global_lock();
          }
          solvep1(kc6[curm], pc6[curm], d, 0, ALLMOVEMASK, CANONSEQSTART);
        }
    }

    moveseq sol;
    int m = cubepos::invm[(solmap % 3) * KOCSYMM];
    for (int i = 0; i < bestsol; i++)
      sol.push_back(cubepos::move_map[m][bestmoves[i]]);
    if (solmap >= 3)
      sol = cubepos::invert_sequence(sol);
    cubepos cpt;
    for (unsigned int i = 0; i < sol.size(); i++)
      cpt.move(sol[i]);
    if (cpt != pos)
      error("! move sequence doesn't work");
    return solution(cp, seq, phase2probes, sol);
  }

  void solvep1(const kocsymm &kc, const permcube &pc, int togo, int sofar,
               int movemask, int canon) {
    if (togo == 0) {
      if (kc == identity_kc)
        solvep2(pc, sofar);
      return;
    }
    if (finished)
      return;
    togo--;
    kocsymm kc2;
    permcube pc2;
    int newmovemask;
    while (!finished && movemask) {
      int mv = ffs1(movemask);
      movemask &= movemask - 1;
      kc2 = kc;
      kc2.move(mv);
      int nd = phase1prune::lookup(kc2, togo, newmovemask);
      if (nd <= togo && (togo == nd || togo + nd >= 5)) {
        pc2 = pc;
        pc2.move(mv);
        moves[sofar] = mv;
        int new_canon = cubepos::next_cs(canon, mv);
        solvep1(kc2, pc2, togo, sofar + 1,
                newmovemask & cubepos::cs_mask(new_canon), new_canon);
      }
    }
  }

  void solvep2(const permcube &pc, int sofar) {
    phase2probes++;
    int d = phase2prune::lookup(pc);
    if (d + sofar <= keepbound) {
      moveseq ms = phase2prune::solve(pc, keepbound - sofar);
      if ((int)(ms.size()) + sofar <= keepbound &&
          (ms.size() > 0 || pc == identity_pc)) {
        int cursol = ms.size() + sofar;
        for (unsigned int i = 0; i < ms.size(); i++)
          moves[sofar + i] = ms[i];
        if (cursol < bestsol) {
          bestsol = cursol;
          memcpy(bestmoves, moves, bestsol);
          if (verbose > 1) {
            get_global_lock();
            cout << "New solution for " << seq << " at " << bestsol << endl;
            release_global_lock();
          }
          solmap = curm;
        }
        if (numsols > 1) {
          get_global_lock();
          moveseq sol;
          int m = cubepos::invm[(curm % 3) * KOCSYMM];
          for (int i = 0; i < cursol; i++)
            sol.push_back(cubepos::move_map[m][moves[i]]);
          if (curm >= 3)
            sol = cubepos::invert_sequence(sol);
          cout << "TSOL " << cursol << " " << cubepos::moveseq_string(sol)
               << endl;
          release_global_lock();
          keepcounts[cursol]++;
          keepsum++;
          while (keepbound > 0 && keepsum >= numsols) {
            keepsum -= keepcounts[keepbound];
            keepbound--;
          }
        } else {
          keepbound = bestsol - 1;
        }
        if (bestsol <= target_length)
          finished = 1;
      }
    }
    if (phase2probes >= phase2limit && bestsol < MAX_MOVES)
      finished = 1;
  }
};
