#include "2sidekmer.h"

static int strict_contains_neighbours( char * query, struct bloom *sparse_bloom,
      uint8_t dist, uint8_t left, size_t kmer_size, uint8_t s, char * neighbour,
      int * final_contain );

int parse_fasta ( FILE * fd, struct bloom * bloom, struct bloom * edge_bloom,
        size_t kmer_size )
{
  char kmer[kmer_size + 1];
  char *line = NULL, last_kmer[kmer_size + 1], edge_kmer[kmer_size + 1];
  size_t len = 0;
  ssize_t read;
  uint8_t left_edge = 0, first = 1, i;
  uint8_t check;
  unsigned int line_cnt = 0;

  // printf("\n\n******** Parse FASTA ********\n\n");

  memset(kmer, 0, sizeof(kmer));

  while((read = getline(&line, &len, fd)) != -1)
  {
    if(strstr(line, ">"))
    {
      if(line_cnt)
      {
        bloom_add(edge_bloom, edge_kmer, kmer_size);
        // printf("right edge kmer to add: %s\n", edge_kmer);
      }
      left_edge = 1;
      line_cnt++;
      first = 1;
      continue;
    }
    for (i = 0; i < len - 1; i++)
    {
      if(((strlen(line + i) - 1 < kmer_size)
          || (strlen(line + i) + strlen(last_kmer) - 1 < kmer_size))
          && first)
      {
        // printf("break\n");
        snprintf(last_kmer, kmer_size + 1, "%s", line + i - 1);
        first = 0;
        break;
      }
      memset(kmer, 0, sizeof(kmer));
      if(first)
      {
        snprintf(kmer, kmer_size+1, "%s", line+i);
        // printf("in first kmer to add: %s\n", kmer);
        check = bloom_add(bloom, kmer, kmer_size);
        snprintf(last_kmer, kmer_size + 1, "%s", line + i - 1);
        snprintf(edge_kmer, kmer_size + 1, "%s", kmer);
      }
      else
      {
        unsigned int j = 0;
        // printf("in else, last_kmer = %s\n", last_kmer);
        if(strlen(line) < kmer_size)
        {
          /* If length of line is < kmer_size (it could be the last line of
             sequence read) don't include new line character '\n' */
          *(line + strlen(line) - 1) = '\0';
        }
        while(strlen(last_kmer) > 0)
        {
          memset(kmer, 0, sizeof(kmer));
          j++;
          memmove(last_kmer, &(last_kmer[1]), strlen(&(last_kmer[1])));
          last_kmer[kmer_size - j] = '\0';
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
          check = bloom_add(bloom, kmer, kmer_size);
          snprintf(edge_kmer, kmer_size + 1, "%s", kmer);
        }
        first = 1;
      }
      if(left_edge)
      {
        // printf("left edge kmer to add: %s\n", kmer);
        check = bloom_add(edge_bloom, kmer, kmer_size);
        left_edge = 0;
      }
    }
    first = 0;
  }
  // printf("right edge kmer to add: %s\n", kmer);
  check = bloom_add(edge_bloom, kmer, kmer_size);
  if(line)
  {
    free(line);
  }
  // printf("\n\n***********************************\n\n");
  return 0;
}

int two_sided_contains ( char * kmer, struct bloom * bloom,
      struct bloom * edge_bloom, size_t kmer_size )
{
  char bases[] = {'A', 'C', 'G', 'T'};
  char tmp_kmer[kmer_size+1];
  uint8_t i, contains_left = 0, contains_right = 0, check;
  // printf("******* Two-sided contains ********\n");
  for (i = 0; i < sizeof(bases); i++)
  {
    memset(tmp_kmer, 0, sizeof(tmp_kmer));
    snprintf(tmp_kmer, strlen(kmer), "%s", kmer+1);
    strncat(tmp_kmer, &bases[i], 1);
    check = bloom_check(bloom, tmp_kmer, kmer_size);
    if(check)
    {
      // printf("right neighbour %s of %s kmer present\n", tmp_kmer, kmer);
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
    check = bloom_check(bloom, tmp_kmer, kmer_size);
    if(check)
    {
      // printf("left neighbour %s of %s kmer present\n", tmp_kmer, kmer);
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
    if(bloom_check(edge_bloom, kmer, kmer_size))
    {
      // printf("kmer %s contained in edges\n", kmer);
      return 1;
    }
  }
  return 0;
}

int sparse_fasta ( FILE *fd, struct bloom * sparse_bloom,
      struct bloom * edge_bloom, size_t kmer_size, uint8_t s )
{
  char kmer[kmer_size + 1];
  char *line = NULL, last_kmer[kmer_size + 1], edge_kmer[kmer_size + 1];
  size_t len = 0;
  ssize_t read;
  uint8_t left_edge = 0, first = 1, i;
  uint8_t check;
  uint8_t cnt = 0;
  unsigned int line_cnt = 0;

  // printf("\n\n******** Sparse FASTA ********\n\n");

  memset(kmer, 0, sizeof(kmer));

  while((read = getline(&line, &len, fd)) != -1)
  {
    if(strstr(line, ">"))
    {
      if(line_cnt)
      {
        bloom_add(edge_bloom, edge_kmer, kmer_size);
        first = 1;
      }
      left_edge = 1;
      line_cnt++;
      continue;
    }
    for (i = 0; i < len - 1; i++)
    {
      if(((strlen(line + i) - 1 < kmer_size)
          || (strlen(line + i) + strlen(last_kmer) - 1 < kmer_size))
          && first)
      {
        snprintf(last_kmer, kmer_size + 1, "%s", line + i - 1);
        first = 0;
        break;
      }
      memset(kmer, 0, sizeof(kmer));
      if(first)
      {
        snprintf(kmer, kmer_size+1, "%s", line+i);
        if(!(cnt % s))
        {
          check = bloom_add(sparse_bloom, kmer, kmer_size);
          // printf("in first kmer to add: %s\n", kmer);
        }
        cnt++;
        if(cnt == s+1)
        {
          cnt = 1;
        }
        snprintf(last_kmer, kmer_size + 1, "%s", line + i - 1);
        snprintf(edge_kmer, kmer_size + 1, "%s", kmer);
      }
      else
      {
        unsigned int j = 0;
        // printf("in else, last_kmer = %s\n", last_kmer);
        if(strlen(line) < kmer_size)
        {
          /* If length of line is < kmer_size (it could be the last line of
             sequence read) don't include new line character '\n' */
          *(line + strlen(line) - 1) = '\0';
        }
        while(strlen(last_kmer) > 0)
        {
          memset(kmer, 0, sizeof(kmer));
          j++;
          memmove(last_kmer, &(last_kmer[1]), strlen(&(last_kmer[1])));
          last_kmer[strlen(last_kmer) - 1] = '\0';
          snprintf(kmer, strlen(last_kmer) + 1, "%s", last_kmer);
          strncat(kmer, line, kmer_size - strlen(last_kmer));
          /* in case that last line of input file has less characters than
            kmer_size, break at last k-mer. If break statement is missing then
            the last k-mer won't be of kmer_size length */
          if(strlen(kmer) < kmer_size)
          {
            break;
          }
          if(!(cnt % s))
          {
            // printf("last_kmer %s\n", last_kmer);
            // printf("in else kmer to add: %s\n", kmer);
            check = bloom_add(sparse_bloom, kmer, kmer_size);
          }
          cnt++;
          if(cnt == s+1)
          {
            cnt = 1;
          }
          snprintf(edge_kmer, kmer_size + 1, "%s", kmer);
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
  check = bloom_add(sparse_bloom, kmer, kmer_size);
  if(line)
  {
    free(line);
  }
  // printf("\n\n***********************************\n\n");
  return 0;
}

static int strict_contains_neighbours ( char * query,
      struct bloom * sparse_bloom, uint8_t dist, uint8_t left, size_t kmer_size,
      uint8_t s, char * neighbour, int * final_contain )
{
  char bases[] = {'A', 'C', 'G', 'T'};
  int contains;
  if( dist == 0)
  {
    if(*final_contain)
    {
      return 1;
    }
    else
    {
      return 0;
    }
  }
  for(int i = 0; i < sizeof(bases); i++)
  {
    if(left)
    {
      neighbour[s-dist] = bases[i];
    }
    else
    {
      neighbour[kmer_size - dist] = bases[i];
    }
    // snprintf(&neighbour[s], kmer_size - s + 1, "%s", query);
    // printf("reconstructed neighbour of %s is %s, with %d dist\n", query,
    //         neighbour, dist);
    if(dist == 1){
      contains = bloom_check(sparse_bloom, neighbour, kmer_size);
      if (contains)
      {
        if(left)
        {
          // printf("left neighbour %s of %s is present in bloom filter\n",
          //         neighbour, query);
        }
        else
        {
          // printf("right neighbour %s of %s is present in bloom filter\n",
          //         neighbour, query);
        }
      }
      *final_contain |= contains;
      if(contains)
      {
        return 1;
      }
    }
    strict_contains_neighbours(query, sparse_bloom, dist - 1, left, kmer_size,
          s, neighbour, final_contain);
  }
  return *final_contain;
}

int strict_contains ( char * query, struct bloom * sparse_bloom,
      struct bloom * edge_bloom, uint8_t s, size_t kmer_size )
{
  char neighbour[kmer_size + 1];
  int contains_left = 0, contains_right = 0, contains = 0;
  uint8_t i;
  if(bloom_check(sparse_bloom, query, kmer_size))
  {
    return 1;
  }
  snprintf(&neighbour[s], kmer_size - s + 1, "%s", query);
  contains_left = strict_contains_neighbours(query, sparse_bloom, s, 1,
        kmer_size, s, neighbour, &contains);
  snprintf(neighbour, kmer_size - s + 1, "%s", query+s);
  contains = 0;
  contains_right = strict_contains_neighbours(query, sparse_bloom, s, 0,
        kmer_size, s, neighbour, &contains);
  if (contains_left && contains_right)
  {
    // printf("contains and\n");
    return 1;
  }
  if (contains_left || contains_right)
  {
    // printf("contains or\n");
    if(bloom_check(edge_bloom, query, kmer_size))
    {
      return 1;
    }
    if(contains_right)
    {
      // printf("contains right\n");
      if(bloom_check(edge_bloom, neighbour, kmer_size))
      {
        return 1;
      }
    }
  }
  for (i = 0; i < s; i++)
  {
    // printf("********* i = %d ***********\n", i);
    contains = 0;
    snprintf(&neighbour[i], kmer_size - i + 1, "%s", query);
    // printf("first naighbour %s of %s for %d\n", neighbour, query, i);
    contains_left = strict_contains_neighbours(query, sparse_bloom, i, 1,
          kmer_size, i, neighbour, &contains);
    // printf("****************************\n");
    snprintf(neighbour, kmer_size - (s - (i + 1)), "%s", query+(s - (i)));
    // printf("first naighbour %s of %s for %d\n", neighbour, query, i);
    contains = 0;
    contains_right = strict_contains_neighbours(query, sparse_bloom, s - (i), 0,
        kmer_size, s - (i), neighbour, &contains);
    if (contains_left && contains_right)
    {
      return 1;
    }
    else
    {
      if (contains_left || contains_right)
      {
        if(bloom_check(edge_bloom, query, kmer_size))
        {
          return 1;
        }
      }
      if(contains_right)
      {
        // printf("contains right\n");
        if(bloom_check(edge_bloom, neighbour, kmer_size))
        {
          return 1;
        }
      }
    }
  }
  return 0;
}
