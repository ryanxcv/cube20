#pragma once

#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>

static bool read_table(unsigned char **p, const char *filename, int memsize) {
  int fd = open(filename, O_RDONLY);
  if (fd == -1)
    return false;
  *p = (unsigned char *) mmap(NULL, memsize, PROT_READ, MAP_SHARED, fd, 0);
  if (*p == MAP_FAILED) {
      close(fd);
      cerr << "mmap failed" << endl;
      return false;
  }
  return true;
}

static inline void read_table_req(unsigned char **p, const char *filename, int memsize) {
  if (not read_table(p, filename, memsize)) {
    cerr << "Error: failed to read table: " << filename << endl;
    exit(-1);
  }
}

static inline void write_table(const char *filename, unsigned char *mem, unsigned int memsize, int file_checksum = 0) {
  FILE *f = fopen(filename, "wb");
  if (f == 0)
    error("! cannot write pruning file to current directory");
  if (fwrite(mem, 1, memsize, f) != memsize)
    error("! error writing pruning table");
  if (fwrite(&file_checksum, sizeof(int), 1, f) != 1)
    error("! error writing pruning table");
  fclose(f);
}
