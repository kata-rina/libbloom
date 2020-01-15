#include "2sidekmer.h"
#include <bloom.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <mutate.h>
#include <sys/time.h>
#include <time.h>
#include <best_fit.h>
#include<hitting_set.h>
#include<relaxed_sparse.h>
#include<onesidekbf.h>

#define S_DIST          1

/*
    ~ Constant K indicates length of k-mers.
    !!! needs to be user defined: TO DO !!!
*/
#define K       20

int million_two_sided ( FILE * fd, size_t kmer_size );
int million_strict ( FILE * fd, size_t kmer_size, uint8_t s );
int million_bloom ( FILE * fd, size_t kmer_size );

int test_kbf(struct bloom * bloom_filter, FILE *f, int kmer_size, int type);
float calculate_fpr(int queries, int positive, int mutated);
int test_relaxed(struct bloom *bloom, struct bloom *edge, int kmer_size, FILE *f, int s);

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

unsigned int queries = QUERIES, present = 0, strict_present = 0;
unsigned int bloom_present = 0;

/* ~ main function can have up to 4 input arguments:
          1. path to file to test upon
          2. size of one k-mer (k)
          3. s distance
          4. 0 - don't mutate, 1 - mutate
  ~ if there isn't enough arguments, code will process data with hard coded
    constants, k will be 20 and s will be 1.
  ~ usage example: ./test ../test_files/test.txt 20 1 1 */

int main (int argc, char ** argv)
{
  FILE * fd, *query_fd;
  uint64_t fsize;
  unsigned int mutated = 0;
  int dist = S_DIST;
  int k = K, mut = 1;
  // char *kmer = "AAGAGACCGGCGATTCTAGT";
  long init_before, init_after, parse_before, parse_after;
  long query_before, query_after;
  int parse_time;
  struct timeval tp;
  if(argc > 1)
  {
    fd = fopen(argv[1], "r");
  }
  else{
    fd = fopen("../test_files/genomic.txt", "r");
  }
  if(!fd)
  {
    printf("Unable to open file\n");
    return 0;
  }
  if(argc >=3)
  {
    k = atoi(argv[2]);
    if(argc > 3)
    {
      dist = atoi(argv[3]);
    }
    if(argc > 4)
    {
      mut = atoi(argv[4]);
      if (mut != 1 && mut != 0)
      {
        printf("Argument 4 must be 0 or 1\n");
        return 0;
      }
    }
  }
  mutated = mutate(fd, k, mut);
  query_fd = fopen("../test_files/mutate.txt", "r");
  if(!query_fd)
  {
    printf("Cannot open file\n");
    return 0;
  }
  fseek(fd, 0, SEEK_END);
  fsize = ftell(fd);
  if(fsize < 1000)
  {
    fsize *= 100000;
  }
  fsize *= 4;
  // printf("Size is %ld\n", fsize);
  gettimeofday(&tp, NULL);
  long before = (tp.tv_sec * 1000L) + (tp.tv_usec / 1000L);
  fseek(fd, 0, SEEK_SET);
  init_before = (tp.tv_sec * 1000L) + (tp.tv_usec / 1000L);
  if(bloom_init(&bloom, fsize, 0.29))
  {
    printf("Bloom not initialized\n");
    return 0;
  }
  if(bloom_init(&edge_bloom, fsize, 0.28))
  {
    printf("Edge bloom not initialized\n");
    return 0;
  }
  printf("Size of bloom filter is %ld\n", fsize);
  gettimeofday(&tp, NULL);
  init_after = (tp.tv_sec * 1000L) + (tp.tv_usec / 1000L);
  printf("********************Simple bloom query*********************\n");
  parse_before = init_after;
  parse_fasta(fd, &bloom, &edge_bloom, k);
  gettimeofday(&tp, NULL);
  parse_after = (tp.tv_sec * 1000L) + (tp.tv_usec / 1000L);
  fseek(query_fd, 0, SEEK_SET);
  gettimeofday(&tp, NULL);
  query_before = (tp.tv_sec * 1000L) + (tp.tv_usec / 1000L);
  million_bloom(query_fd, k);
  gettimeofday(&tp, NULL);
  query_after = (tp.tv_sec * 1000L) + (tp.tv_usec / 1000L);
  gettimeofday(&tp, NULL);
  long after = (tp.tv_sec * 1000L) + (tp.tv_usec / 1000L);
  parse_time = (int)(init_after - init_before + parse_after - parse_before);
  printf("  Building bloom in %d ms\n  Query time is %d ms\n", parse_time,
      (int)(query_after - query_before));
  printf("FPR is %f\n", (float)(-(queries-mutated-bloom_present)/(float)queries));
  printf("%d k-mers present of %d standard queries with %d mutated\n",
      bloom_present, queries, mutated);
  printf("**********************2-sided query***********************\n");
  fseek(query_fd, 0, SEEK_SET);
  gettimeofday(&tp, NULL);
  before = (tp.tv_sec * 1000L) + (tp.tv_usec / 1000L);
  million_two_sided(query_fd, k);
  gettimeofday(&tp, NULL);
  after = (tp.tv_sec * 1000L) + (tp.tv_usec / 1000L);
  // parse_time = (int) (init_after - init_before + parse_after - parse_before +
  //     after - before);
  printf("  Building bloom in %d ms\n  Query time is %d ms\n", parse_time,
      (int) (after - before));
  printf("FPR is %f\n", (float)(-(queries-mutated-present)/(float)queries));
  printf("%d k-mers present of %d queries in bloom with %d mutated\n", present,
      queries, mutated);
  printf("***********************Strict query************************\n");
  gettimeofday(&tp, NULL);
  init_before = (tp.tv_sec * 1000L) + (tp.tv_usec / 1000L);
  if(bloom_init(&sparse_bloom, fsize/1.5, 0.28))
  {
    printf("Sparse bloom not initialized\n");
    return 0;
  }
  gettimeofday(&tp, NULL);
  init_after = (tp.tv_sec * 1000L) + (tp.tv_usec / 1000L);
  fseek(fd, 0, SEEK_SET);
  gettimeofday(&tp, NULL);
  parse_before = (tp.tv_sec * 1000L) + (tp.tv_usec / 1000L);
  // sparse_fasta(fd, &sparse_bloom, &edge_bloom, k, dist);
  best_fit_parse(fd, &sparse_bloom, &edge_bloom, k, dist);
  gettimeofday(&tp, NULL);
  parse_after = (tp.tv_sec * 1000L) + (tp.tv_usec / 1000L);
  fseek(query_fd, 0, SEEK_SET);
  gettimeofday(&tp, NULL);
  before = (tp.tv_sec * 1000L) + (tp.tv_usec / 1000L);
  million_strict(query_fd, k, dist);
  gettimeofday(&tp, NULL);
  after = (tp.tv_sec * 1000L) + (tp.tv_usec / 1000L);
  parse_time = (int)(init_after - init_before + parse_after - parse_before);
  printf("  Building bloom in %d ms\n  Query time is %d ms\n", parse_time,
      (int)(after - before));
  printf("FPR is %f\n", (float)(-(queries-mutated-strict_present)/(float)queries));
  printf("%d k-mers present of %d strict queries with %d mutated\n",
      strict_present, queries, mutated);


  // ***************************************************************************
  struct timeval tstart, tstop;
  printf("********** One sided kbf performance testing with mutation **********\n" );
  fseek(query_fd, 0, SEEK_SET);
  gettimeofday(&tstart, NULL);
  int match_mut = test_kbf(&bloom, query_fd , k, 1);
  gettimeofday(&tstop, NULL);

  printf("* Number of mutated k-mers => %d\n", mutated );
  printf("* Positive queries => %d \n", match_mut );
  printf("* Operating time = %f s\n ",
            (double) (tstop.tv_usec - tstart.tv_usec) / 1000000 +
            (double) (tstop.tv_sec - tstart.tv_sec));

  float fpr = calculate_fpr(QUERIES, match_mut, mutated);
  printf("* False positive rate => %f\n", fpr);
  fseek(query_fd, 0, SEEK_SET);
  struct bloom kata_bloom;
  printf("********** Hitting set testing with mutation **********\n" );
  bloom_init(&kata_bloom, fsize, 0.29);

  struct bloom kata_edge_bloom;
  bloom_init(&kata_edge_bloom, 20000, 0.29);

  fseek(fd, 0, SEEK_SET);

  gettimeofday(&tstart, NULL);
  parse_hitting_set(k, fd, &kata_bloom, &kata_edge_bloom);
  gettimeofday(&tstop, NULL);


  printf("* Operating time = %f s\n ",
            (double) (tstop.tv_usec - tstart.tv_usec) / 1000000 +
            (double) (tstop.tv_sec - tstart.tv_sec));

  fseek(query_fd, 0, SEEK_SET);
  // fseek(f_nonmutated, 0, SEEK_SET);
  // int s = 1;
  // int kmer_size = 20;

  gettimeofday(&tstart, NULL);
  int a = test_relaxed( &kata_bloom, &kata_edge_bloom, k, query_fd, dist );
  gettimeofday(&tstop, NULL);
  printf("* Operating time = %f s\n ",
            (double) (tstop.tv_usec - tstart.tv_usec) / 1000000 +
            (double) (tstop.tv_sec - tstart.tv_sec));

  fpr = calculate_fpr(QUERIES, a, mutated);

  printf("False positive rate => %f\n", fpr );
  // ***************************************************************************
  fclose(fd);
  if(query_fd)
  {
    fclose(query_fd);
  }
  fclose(query_fd);
  bloom_free(&kata_bloom);
  bloom_free(&kata_edge_bloom);
  bloom_free(&bloom);
  bloom_free(&edge_bloom);
  bloom_free(&sparse_bloom);
  return 0;
}

// function for testing relaxed_sparse kbf performance
int test_relaxed(struct bloom *bloom, struct bloom *edge, int kmer_size, FILE *f, int s){


  ssize_t read;
  size_t len = 0;
  char *line = NULL;
  char sequence[kmer_size];
  int found = 0; // counter of positive queries
  int cnt = 0;

  // bloom_print(edge);
  // bloom_print(bloom);



  while( ( read = getline( &line, &len, f ) ) != -1 ){

    *(line + kmer_size) = '\0';

    if (relaxed_contains( line, kmer_size, s, bloom, edge ))
    {
      found++;
    }
    else{
      // printf("%s\n", line );
    }
  }

  free( line );
  return found;

}



// function for testing kmer bloom filter performance
// if type => 1 - test one sided KBF
// if type => 0 - test regular KBF
int test_kbf(struct bloom * bloom_filter, FILE *f, int kmer_size, int type){

  ssize_t read;
  size_t len = 0;
  char *line = NULL;
  char sequence[kmer_size];
  int found = 0; // counter of positive queries

  memset( sequence, 0, kmer_size*sizeof(char) );

  while( ( read = getline( &line, &len, f ) ) != -1 ){

    snprintf( sequence, kmer_size + 1, "%s", line );

    if(type == 1){
      if (onesided_kbf_contains( sequence, kmer_size, bloom_filter )){found++;}
    }
    else{
      if (bloom_check(bloom_filter, sequence, kmer_size ) == 1){found++;}
    }

  }

  free( line );
  return found;

}

// returns false positive rate in percantage
float calculate_fpr(int queries, int positive, int mutated){

  float fpr;
  fpr = (float)(((positive + mutated) - queries) / (float)queries);
  // fpr *= 100;
  return fpr;
}

int million_two_sided ( FILE * fd, size_t kmer_size )
{
  // char kmer[kmer_size + 1];
  char *line = NULL;
  size_t len = 0;
  while(getline(&line, &len, fd) != -1)
  {
    *(line + kmer_size) = '\0';
    if(two_sided_contains(line, &bloom, &edge_bloom, kmer_size))
    {
      present++;
    }
    // else
    // {
    //   printf("%s kmer not present\n", line);
    // }
  }
  free(line);
  return present;
}

int million_strict ( FILE *fd, size_t kmer_size, uint8_t s )
{
  char *line = NULL;
  size_t len = 0;
  while(getline(&line, &len, fd) != -1)
  {
    *(line + kmer_size) = '\0';
    if(strict_contains(line, &sparse_bloom, &edge_bloom, s, kmer_size))
    {
      strict_present++;
    }
    // else
    // {
    //   printf("%s kmer not present\n", line);
    // }
  }
  free(line);
  return strict_present;
}

int million_bloom ( FILE *fd, size_t kmer_size )
{
  char *line = NULL;
  size_t len = 0;
  while(getline(&line, &len, fd) != -1)
  {
    *(line + kmer_size) = '\0';
    if(bloom_check(&bloom, line, kmer_size))
    {
      bloom_present++;
    }
    // else
    // {
    //   printf("%s kmer not present\n", line);
    // }
  }
  free(line);
  return bloom_present;
}
