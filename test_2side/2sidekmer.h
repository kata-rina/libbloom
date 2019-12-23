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
    ~ Constant K indicates length of k-mers.
    !!! needs to be user defined: TO DO !!!
*/
#define K       20

/*
    ~ Structure bloom is used for basic bloom-filter which contains all k-mers
    contained in genome sequences.
    ~ Structure edge_bloom is used for bloom-filter which contains only edge
    k-mers of genome sequences.
    ~ Structure sparse_bloom is used to store k-mers with distance of s
*/
struct bloom bloom;
struct bloom edge_bloom;
struct bloom sparse_bloom;

/*
    Funtion parse_fasta is used to decompose input FASTA file in k-mers with
    length of kmer_size and store them in bloom-filters defined above.

    *inputs:  fd - file descriptor of FASTA input
              kmer_size - length of k-mers to be stored in bloom-filters

    It supports only whole genome and not sequence reads:
            !!! TO DO !!!
*/
int parse_fasta ( FILE * fd, size_t kmer_size );
int two_sided_contains ( char * kmer, size_t kmer_size );
int sparse_fasta ( FILE *fd, size_t kmer_size, uint8_t s );
int strict_contains ( char * query, uint8_t s, size_t kmer_size );

#endif
