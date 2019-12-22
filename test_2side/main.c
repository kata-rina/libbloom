#include "2sidekmer.h"
#include <bloom.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

int main (void)
{
  FILE * fd;
  uint64_t fsize;
  unsigned int check = 0;
  fd = fopen("../test_files/test.txt", "r");
  if(!fd)
  {
    printf("Unable to open file\n");
    return 0;
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
  if(bloom_init(&edge_bloom, 200*K, 0.0001))
  {
    printf("Edge bloom not initialized\n");
    return 0;
  }
  if(bloom_init(&sparse_bloom, fsize/10, 0.2))
  {
    printf("Sparse bloom not initialized\n");
    return 0;
  }
  parse_fasta(fd, K);
  fseek(fd, 0, SEEK_SET);
  sparse_fasta(fd, K, 1);
  bloom_print(&bloom);
  bloom_print(&edge_bloom);
  bloom_print(&sparse_bloom);
  // check = bloom_check(&bloom, "TTTTTTTATATATAGGGGCC", K);
  // if(check == 1)
  // {
  //   printf("TTTTTTTATATATAGGGGCC present in bloom\n");
  // }
  // check = bloom_check(&edge_bloom, "TTTTTTTATATATAGGGGCC", K);
  // if(check == 1)
  // {
  //   printf("TTTTTTTATATATAGGGGCC present in edge bloom\n");
  // }
  // check = two_sided_contains("TTGAGGTCGCAGTGACCCCG", K);
  // if(check == 1)
  // {
  //   printf("TTGAGGTCGCAGTGACCCCG contained in bloom filter\n");
  // }
  // check = two_sided_contains("TCATGATTCGGTACCTGGGT", K);
  // if(check == 1)
  // {
  //   printf("TCATGATTCGGTACCTGGGT contained in bloom filter\n");
  // }
  // check = two_sided_contains("TTAAAGAGACCGGCGATTCT", K);
  // if(check == 1)
  // {
  //   printf("TTAAAGAGACCGGCGATTCT contained in bloom filter\n");
  // }
  char neighbour[K+1];
  int dist = 1;
  snprintf(&neighbour[dist], K - dist +1, "%s", "TTAAAGAGACCGGCGATTCT");
  check = strict_contains_neighbours ("TTAAAGAGACCGGCGATTCT", dist, 1, K, dist, neighbour, 0);
  if(check)
  {
    printf("TTAAAGAGACCGGCGATTCT present in sparse bloom\n");
  }
  fclose(fd);
  bloom_free(&bloom);
  bloom_free(&edge_bloom);
  bloom_free(&sparse_bloom);
  return 0;
}
