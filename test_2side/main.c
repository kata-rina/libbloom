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
  fd = fopen("../test_files/genomic.txt", "r");
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
  if(bloom_init(&edge_bloom, 200*K, 0.01))
  {
    printf("Edge bloom not initialized\n");
    return 0;
  }
  parse_fasta(fd, K);
  bloom_print(&bloom);
  bloom_print(&edge_bloom);
  check = bloom_check(&bloom, "TTTTTTTATATATAGGGGCC", K);
  if(check == 1)
  {
    printf("TTTTTTTATATATAGGGGCC present in bloom\n");
  }
  check = bloom_check(&edge_bloom, "TTTTTTTATATATAGGGGCC", K);
  if(check == 1)
  {
    printf("TTTTTTTATATATAGGGGCC present in edge bloom\n");
  }
  check = two_sided_contains("CTGCTTTTATTAAGGTCTTG", K);
  if(check == 1)
  {
    printf("CTGCTTTTATTAAGGTCTTG contained in bloom filter\n");
  }
  fclose(fd);
  bloom_free(&bloom);
  bloom_free(&edge_bloom);
  return 0;
}
