#ifndef RELAXED_SPARSE_H
#define RELAXED_SPARSE_H

#include<bloom.h>
#include<stdint.h>

int decide_present(char *query, int size, int contains_l, int contains_r,
                    struct bloom *edge_bloom);

int relaxed_contains(char *query, int size, int s,
                      struct bloom *bloom, struct bloom *egde_bloom );

int relaxed_contains_neighbours(char *query, int left_dist, int right_dist,
            struct bloom * bloom, struct bloom *edge_bloom, int kmer_size);

void left_neighbours_relaxed ( char * query,
      struct bloom * sparse_bloom, int dist, int kmer_size,
      int s, int *contains);

void right_neighbours_relaxed ( char * query,
      struct bloom * sparse_bloom, int dist, int kmer_size,
      int s, int *contains);

#endif
