#include "2sidekmer.h"

int parse_fasta ( FILE * fd, size_t kmer_size )
{
  char kmer[kmer_size];
  char *line = NULL, last_kmer[kmer_size];
  size_t len = 0;
  ssize_t read;
  uint8_t left_edge = 0, first = 1, i;
  uint8_t check;

  // printf("\n\n******** Parse FASTA ********\n\n");

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
      if(((strlen(line + i) - 1 < kmer_size)
          || (strlen(line + i) + strlen(last_kmer) - 1 < kmer_size))
          && first)
      {
        // printf("break\n");
        snprintf(last_kmer, kmer_size, "%s", line + i);
        first = 0;
        break;
      }
      memset(kmer, 0, sizeof(kmer));
      if(first)
      {
        snprintf(kmer, kmer_size+1, "%s", line+i);
        // printf("in first kmer to add: %s\n", kmer);
        check = bloom_add(&bloom, kmer, kmer_size);
        snprintf(last_kmer, kmer_size, "%s", line + i);
      }
      else
      {
        unsigned int j = 0;
        // printf("in else, last_kmer = %s\n", last_kmer);
        while(strlen(last_kmer) > 0)
        {
          memset(kmer, 0, sizeof(kmer));
          snprintf(kmer, strlen(last_kmer) + 1, "%s", last_kmer);
          strncat(kmer, line, kmer_size - strlen(last_kmer));
          /* in case that last line of input file has less characters than
            kmer_size, break at last k-mer. If break statement is missing then
            the last k-mer won't be of kmer_size length */
          if(strlen(kmer) < kmer_size)
          {
            break;
          }
          // printf("in else kmer to add: %s\n", kmer);
          check = bloom_add(&bloom, kmer, kmer_size);
          j++;
          memmove(last_kmer, &(last_kmer[1]), strlen(&(last_kmer[1])));
          last_kmer[kmer_size - j - 1] = '\0';
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
  // printf("\n\n***********************************\n\n");
  return 0;
}

int two_sided_contains ( char * kmer, size_t kmer_size )
{
  char bases[] = {'A', 'C', 'G', 'T'};
  char tmp_kmer[kmer_size+1];
  uint8_t i, contains_left = 0, contains_right = 0, check;
  printf("******* Two-sided contains ********\n");
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
    // else
    // {
    //   printf("right neighbour %s of %s kmer not present\n", tmp_kmer, kmer);
    // }
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
    // else
    // {
    //   printf("left neighbour %s of %s kmer not present\n", tmp_kmer, kmer);
    // }
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

int sparse_fasta ( FILE *fd, size_t kmer_size, uint8_t s )
{
  char kmer[kmer_size];
  char *line = NULL, last_kmer[kmer_size];
  size_t len = 0;
  ssize_t read;
  uint8_t left_edge = 0, first = 1, i;
  uint8_t check;
  uint8_t cnt = 0;

  // printf("\n\n******** Sparse FASTA ********\n\n");

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
      if(((strlen(line + i) - 1 < kmer_size)
          || (strlen(line + i) + strlen(last_kmer) - 1 < kmer_size))
          && first)
      {
        snprintf(last_kmer, kmer_size, "%s", line + i);
        first = 0;
        break;
      }
      memset(kmer, 0, sizeof(kmer));
      if(first)
      {
        snprintf(kmer, kmer_size+1, "%s", line+i);
        if(!cnt)
        {
          check = bloom_add(&sparse_bloom, kmer, kmer_size);
          // printf("in first kmer to add: %s\n", kmer);
        }
        cnt++;
        if(cnt == s+1)
        {
          cnt = 0;
        }
        snprintf(last_kmer, kmer_size, "%s", line + i);
      }
      else
      {
        unsigned int j = 0;
        // printf("in else, last_kmer = %s\n", last_kmer);
        while(strlen(last_kmer) > 0)
        {
          memset(kmer, 0, sizeof(kmer));
          snprintf(kmer, strlen(last_kmer) + 1, "%s", last_kmer);
          strncat(kmer, line, kmer_size - strlen(last_kmer));
          /* in case that last line of input file has less characters than
            kmer_size, break at last k-mer. If break statement is missing then
            the last k-mer won't be of kmer_size length */
          if(strlen(kmer) < kmer_size)
          {
            break;
          }
          if(!cnt)
          {
            // printf("in else kmer to add: %s\n", kmer);
            check = bloom_add(&sparse_bloom, kmer, kmer_size);
          }
          cnt++;
          if(cnt == s+1)
          {
            cnt = 0;
          }
          j++;
          memmove(last_kmer, &(last_kmer[1]), strlen(&(last_kmer[1])));
          last_kmer[kmer_size - j - 1] = '\0';
        }
        first = 1;
      }
      if(left_edge)
      {
        // printf("left edge kmer to add: %s\n", kmer);
        // check = bloom_add(&edge_bloom, kmer, kmer_size);
        left_edge = 0;
      }
    }
    first = 0;
  }
  // printf("right edge kmer to add: %s\n", kmer);
  // check = bloom_add(&edge_bloom, kmer, kmer_size);
  if(!line)
  {
    free(line);
  }
  // printf("\n\n***********************************\n\n");
  return 0;
}

int strict_contains_neighbours ( char * query, uint8_t dist, uint8_t left,
                                size_t kmer_size, uint8_t s, char * neighbour, int final_contain )
{
  char bases[] = {'A', 'C', 'G', 'T'};
  int contains;
  // int final_contain;
  // char neighbour[kmer_size + 1];
  if( dist == 0 )
  {
    return 0;
  }
  if(left)
  {
    for(int i = 0; i < sizeof(bases); i++)
    {
      // snprintf(&neighbour[s], kmer_size - s + 1, "%s", query);
      neighbour[s-dist] = bases[i];
      printf("reconstructed neighbour of %s is %s, with %d dist\n", query,
              neighbour, dist);
      if(dist == 1){
        contains = bloom_check(&sparse_bloom, neighbour, kmer_size);
        final_contain |= contains;
        // strict_contains_neighbours(query, dist - 1, left,
        //                     kmer_size, s, neighbour, final_contain);
      }
      strict_contains_neighbours(query, dist - 1, left, kmer_size, s, neighbour, final_contain);
    }
    return final_contain;
  }
  else
  {

  }
  return 0;
}
