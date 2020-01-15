#ifndef HITTING_SET_H
#define HITTING_SET_H

#include<bloom.h>
#include<stdio.h>
#include<onesidekbf.h>

typedef struct kmer_node {
  char *previous_kmers;
  char *current_kmer;
  char *next_kmers;

  int previous_count;
  int next_count;
  struct kmer_node *next;

} kmer_node_t;


kmer_node_t *parse_hitting_set(int kmer_size, FILE *f, struct bloom * bloom,
                                struct bloom * edge_bloom);

void add_to_list(char *kmer, char *left, char *right, int kmer_size, kmer_node_t *head);
int check_presence(char *data_set, char *query, int query_length);


#endif
