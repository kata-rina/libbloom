/*
================================================================================
        Header file containing function declarations for two sided k-mer
                      bloom-filter implemetation
================================================================================

    @author: Magdalena Halusek
*/
#ifndef TWOSIDEKMER_H
#define TWOSIDEKMER_H

#define _GNU_SOURCE
#include <bloom.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>

/*
    Funtion parse_fasta is used to decompose input FASTA file in k-mers with
    length of kmer_size and store them in bloom-filters defined above.

    *inputs:  fd - file descriptor of FASTA input
              kmer_size - length of k-mers to be stored in bloom-filters

    It supports only whole genome and not sequence reads:
            !!! TO DO !!!
*/
int parse_fasta ( FILE * fd, struct bloom * bloom, struct bloom * edge_bloom,
      size_t kmer_size );
int two_sided_contains ( char * kmer, struct bloom * bloom,
      struct bloom * edge_bloom, size_t kmer_size );
int sparse_fasta ( FILE *fd, struct bloom * bloom, struct bloom * edge_bloom,
      size_t kmer_size, uint8_t s );
int strict_contains ( char * query, struct bloom * sparse_bloom,
      struct bloom * edge_bloom, uint8_t s, size_t kmer_size );

#endif
