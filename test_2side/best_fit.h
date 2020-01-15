#ifndef BEST_FIT_H
#define BEST_FIT_H

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <bloom.h>
#include <stdlib.h>

int best_fit_parse ( FILE *fd, struct bloom * sparse_bloom,
      struct bloom * edge_bloom, size_t kmer_size, uint8_t s );

#endif
