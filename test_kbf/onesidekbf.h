#ifndef ONESIDEKBF_H
#define ONESIDEKBF_H

#include<stdio.h>
#include<bloom.h>

#define KMER_SIZE   20
#define KBF_SIZE    100000000

// int parse_fasta(FILE * fp, int kmer_size, struct bloom * bloom);
int onesided_kbf_contains( char *kmer, int kmer_size, struct bloom * bloom );
#endif
