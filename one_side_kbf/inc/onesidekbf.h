#ifndef ONESIDEKBF_H
#define ONESIDEKBF_H

#include<stdio.h>
#include<bloom.h>

int genom_add(FILE * fp, int kmer_size, struct bloom * bloom);
int onesided_kbf_check( char *kmer, int kmer_size, struct bloom * bloom );
#endif
