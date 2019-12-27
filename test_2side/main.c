#include "2sidekmer.h"
#include <bloom.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <mutate.h>

#define S_DIST          1

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

unsigned int queries = 0, present = 0, strict_present = 0;

/* ~ main function can have up to 4 input arguments:
          1. path to file to test upon
          2. size of one k-mer (k)
          3. s distance
          4. query to find
  ~ if there isn't enough arguments, code will process data with hard coded
    constants, k will be 20 and s will be 1.
  ~ usage example: ./test ../test_files/test.txt 20 1 "AAGAGACCGGCGATTCTAGT" */

int main (int argc, char ** argv)
{
  FILE * fd, *query_fd;
  uint64_t fsize;
  unsigned int check = 0;
  int dist = S_DIST;
  int k = K;
  char *kmer = "AAGAGACCGGCGATTCTAGT";
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
      // kmer = argv[4];
      query_fd = fopen(argv[4], "r");
    }
  }
  // else
  // {
  //   query_fd = fopen("../test_files/mil_query.txt", "r");
  // }
  // if(!query_fd)
  // {
  //   printf("Unable to open query file\n");
  //   return 0;
  // }
  queries = mutate(fd, k);
  fseek(fd, 0, SEEK_END);
  fsize = ftell(fd);
  if(fsize < 1000)
  {
    fsize *= 1000;
  }
  // fsize = 10000000;
  printf("Size is %ld\n", fsize);
  fseek(fd, 0, SEEK_SET);
  if(bloom_init(&bloom, fsize, 0.000001))
  {
    printf("Bloom not initialized\n");
    return 0;
  }
  if(bloom_init(&edge_bloom, 200*k, 0.000001))
  {
    printf("Edge bloom not initialized\n");
    return 0;
  }
  if(bloom_init(&sparse_bloom, fsize/10, 0.000001))
  {
    printf("Sparse bloom not initialized\n");
    return 0;
  }
  parse_fasta(fd, &bloom, &edge_bloom, k);
  fseek(fd, 0, SEEK_SET);
  sparse_fasta(fd, &sparse_bloom, &edge_bloom, k, dist);
  bloom_print(&bloom);
  bloom_print(&edge_bloom);
  bloom_print(&sparse_bloom);
  // check = bloom_check(&bloom, "TTTTTTTATATATAGGGGCC", k);
  // if(check == 1)
  // {
  //   printf("TTTTTTTATATATAGGGGCC present in bloom\n");
  // }
  // check = bloom_check(&edge_bloom, "TTTTTTTATATATAGGGGCC", k);
  // if(check == 1)
  // {
  //   printf("TTTTTTTATATATAGGGGCC present in edge bloom\n");
  // }
  // check = two_sided_contains("TTGAGGTCGCAGTGACCCCG", &bloom, &edge_bloom, k);
  // if(check == 1)
  // {
  //   printf("TTGAGGTCGCAGTGACCCCG contained in bloom filter\n");
  // }
  // check = two_sided_contains("TCATGATTCGGTACCTGGGT", &bloom, &edge_bloom, k);
  // if(check == 1)
  // {
  //   printf("TCATGATTCGGTACCTGGGT contained in bloom filter\n");
  // }
  // check = two_sided_contains("TTAAAGAGACCGGCGATTCT", &bloom, &edge_bloom, k);
  // if(check == 1)
  // {
  //   printf("TTAAAGAGACCGGCGATTCT contained in bloom filter\n");
  // }
  check = strict_contains (kmer, &sparse_bloom, &edge_bloom, dist, k);
  if(check)
  {
    printf("%s present in sparse bloom\n", kmer);
  }
  fclose(fd);
  // fclose(query_fd);
  bloom_free(&bloom);
  bloom_free(&edge_bloom);
  bloom_free(&sparse_bloom);
  return 0;
}
