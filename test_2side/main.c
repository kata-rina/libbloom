#include "2sidekmer.h"
#include <bloom.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#define S_DIST          1

/*
    ~ Constant K indicates length of k-mers.
    !!! needs to be user defined: TO DO !!!
*/
#define K       20

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
  FILE * fd;
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
      kmer = argv[4];
    }
  }
  fseek(fd, 0, SEEK_END);
  fsize = ftell(fd);
  fsize *= 1000;
  // fsize = 10000000;
  printf("Size is %ld\n", fsize);
  fseek(fd, 0, SEEK_SET);
  if(bloom_init(&bloom, fsize, 0.61))
  {
    printf("Bloom not initialized\n");
    return 0;
  }
  if(bloom_init(&edge_bloom, 200*k, 0.0001))
  {
    printf("Edge bloom not initialized\n");
    return 0;
  }
  if(bloom_init(&sparse_bloom, fsize/10, 0.2))
  {
    printf("Sparse bloom not initialized\n");
    return 0;
  }
  parse_fasta(fd, k);
  fseek(fd, 0, SEEK_SET);
  sparse_fasta(fd, k, dist);
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
  // check = two_sided_contains("TTGAGGTCGCAGTGACCCCG", k);
  // if(check == 1)
  // {
  //   printf("TTGAGGTCGCAGTGACCCCG contained in bloom filter\n");
  // }
  // check = two_sided_contains("TCATGATTCGGTACCTGGGT", k);
  // if(check == 1)
  // {
  //   printf("TCATGATTCGGTACCTGGGT contained in bloom filter\n");
  // }
  // check = two_sided_contains("TTAAAGAGACCGGCGATTCT", k);
  // if(check == 1)
  // {
  //   printf("TTAAAGAGACCGGCGATTCT contained in bloom filter\n");
  // }
  // char neighbour[K+1];
  // snprintf(&neighbour[dist], k - dist +1, "%s", "TTAAAGAGACCGGCGATTCT");
  check = strict_contains (kmer, dist, k);
  if(check)
  {
    printf("%s present in sparse bloom\n", kmer);
  }
  fclose(fd);
  bloom_free(&bloom);
  bloom_free(&edge_bloom);
  bloom_free(&sparse_bloom);
  return 0;
}
