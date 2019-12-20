#include "2sidekmer.h"

int parse_fasta ( FILE * fd, size_t kmer_size )
{
  char kmer[kmer_size];
  char *line = NULL, *last_kmer = NULL;
  size_t len = 0;
  ssize_t read;
  uint8_t left_edge = 0, first = 1, i;
  uint8_t check;

  memset(kmer, 0, sizeof(kmer));

  while((read = getline(&line, &len, fd)) != -1)
  {
    if(strstr(line, ">"))
    {
      left_edge = 1;
      continue;
    }
    for (i = 0; i < len - 1; i++)
    {
      if(strlen(line + i) - 1 < kmer_size)
      {
        last_kmer = line + i + 1;
        break;
      }
      memset(kmer, 0, sizeof(kmer));
      if(first)
      {
        snprintf(kmer, kmer_size+1, "%s", line+i);
        // printf("kmer to add: %s\n", kmer);
        if(strstr(kmer, "GTGCTTTTATTAAGGTCTT"))
        {
          printf("kmer to add: %s\n", kmer);
        }
        check = bloom_add(&bloom, kmer, kmer_size);
        last_kmer = line + i + 1;
      }
      else
      {
        int j = 0;
        while(strlen(last_kmer) > 1)
        {
          memset(kmer, 0, sizeof(kmer));
          snprintf(kmer, strlen(last_kmer) - 1, "%s", last_kmer);
          strncat(kmer, line + j, kmer_size - strlen(last_kmer) + 1);
          if(strstr(kmer, "GTGCTTTTATTAAGGTCTT"))
          {
            printf("kmer to add: %s\n", kmer);
          }
          check = bloom_add(&bloom, kmer, kmer_size);
          last_kmer++;
          j++;
        }
        first = 1;
      }
      if(left_edge)
      {
        // printf("left edge kmer to add: %s\n", kmer);
        check = bloom_add(&edge_bloom, kmer, kmer_size);
        left_edge = 0;
      }
    }
    first = 0;
  }
  // printf("right edge kmer to add: %s\n", kmer);
  check = bloom_add(&edge_bloom, kmer, kmer_size);
  if(!line)
  {
    free(line);
  }
  return 0;
}

int two_sided_contains ( char * kmer, size_t kmer_size )
{
  char bases[] = {'A', 'C', 'G', 'T'};
  char tmp_kmer[kmer_size+1];
  uint8_t i, contains_left = 0, contains_right = 0, check;
  for (i = 0; i < sizeof(bases); i++)
  {
    memset(tmp_kmer, 0, sizeof(tmp_kmer));
    snprintf(tmp_kmer, strlen(kmer), "%s", kmer+1);
    strncat(tmp_kmer, &bases[i], 1);
    check = bloom_check(&bloom, tmp_kmer, kmer_size);
    if(check)
    {
      printf("right neighbour %s of %s kmer present\n", tmp_kmer, kmer);
    }
    else
    {
      printf("right neighbour %s of %s kmer not present\n", tmp_kmer, kmer);
    }
    contains_right |= check;
  }
  for (i = 0; i < sizeof(bases); i++)
  {
    memset(tmp_kmer, 0, sizeof(tmp_kmer));
    tmp_kmer[0] = bases[i];
    strncat(tmp_kmer, kmer, strlen(kmer) - 1);
    check = bloom_check(&bloom, tmp_kmer, kmer_size);
    if(check)
    {
      printf("left neighbour %s of %s kmer present\n", tmp_kmer, kmer);
    }
    else
    {
      printf("left neighbour %s of %s kmer not present\n", tmp_kmer, kmer);
    }
    contains_left |= check;
  }
  if (contains_left && contains_right)
  {
    return 1;
  }
  if(contains_left || contains_right)
  {
    if(bloom_check(&bloom, kmer, kmer_size))
    {
      printf("kmer %s contained in edges\n", kmer);
      return 1;
    }
  }
  return 0;
}
