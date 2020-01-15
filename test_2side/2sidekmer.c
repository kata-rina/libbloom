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
  uint8_t left_edge = 0, first = 1, i, line_length, kmer_length;
  // uint8_t check;
  unsigned int line_cnt = 0, j = 0;

  // printf("\n\n******** Parse FASTA ********\n\n");

  // memset(kmer, 0, sizeof(kmer));
  kmer[kmer_size] = '\0';
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
      line_length = strlen(line + i);
      kmer_length = strlen(last_kmer);
      if(((line_length - 1 < kmer_size)
          || (line_length + kmer_length - 1 < kmer_size))
          && first)
      {
        // printf("break\n");
        snprintf(last_kmer, kmer_size + 1, "%s", line + i - 1);
        first = 0;
        break;
      }
      // memset(kmer, 0, sizeof(kmer));
      if(first)
      {
        snprintf(kmer, kmer_size+1, "%s", line+i);
        // printf("in first kmer to add: %s\n", kmer);
        bloom_add(bloom, kmer, kmer_size);
        // snprintf(last_kmer, kmer_size + 1, "%s", line + i - 1);
        snprintf(edge_kmer, kmer_size + 1, "%s", kmer);
      }
      else
      {
        j = 0;
        line_length = strlen(line);
        // printf("in else, last_kmer = %s\n", last_kmer);
        if(line_length < kmer_size)
        {
          /* If length of line is < kmer_size (it could be the last line of
             sequence read) don't include new line character '\n' */
          *(line + line_length - 1) = '\0';
        }
        kmer_length = strlen(last_kmer);
        while(kmer_length > 0)
        {
          // memset(kmer, 0, sizeof(kmer));
          j++;
          memmove(last_kmer, &(last_kmer[1]), kmer_length - 1);
          last_kmer[kmer_size - j] = '\0';
          kmer_length--;
          snprintf(kmer, kmer_length + 1, "%s", last_kmer);
          kmer[kmer_length] = '\0';
          strncat(kmer, line, kmer_size - kmer_length);
          /* in case that last line of input file has less characters than
            kmer_size, break at last k-mer. If break statement is missing then
            the last k-mer won't be of kmer_size length */
          if(strlen(kmer) < kmer_size)
          {
            break;
          }
          // printf("in else kmer to add: %s\n", kmer);
          bloom_add(bloom, kmer, kmer_size);
          snprintf(edge_kmer, kmer_size + 1, "%s", kmer);
        }
        first = 1;
      }
      if(left_edge)
      {
        // printf("left edge kmer to add: %s\n", kmer);
        bloom_add(edge_bloom, kmer, kmer_size);
        left_edge = 0;
      }
    }
    first = 0;
  }
  // printf("right edge kmer to add: %s\n", kmer);
  bloom_add(edge_bloom, kmer, kmer_size);
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
  char left_kmer[kmer_size+1], right_kmer[kmer_size+1];
  uint8_t i, contains_left = 0, contains_right = 0;
  // printf("******* Two-sided contains ********\n");
  if(!bloom_check(bloom, kmer, kmer_size))
  {
    return 0;
  }
  // memset(tmp_kmer, 0, sizeof(tmp_kmer));
  snprintf(right_kmer, kmer_size, "%s", kmer+1);
  snprintf(&left_kmer[1], kmer_size, "%s" ,kmer);
  left_kmer[kmer_size] = '\0';
  right_kmer[kmer_size] = '\0';
  for (i = 0; i < 4; i++)
  {
    // if(!contains_right)
    // {
      right_kmer[kmer_size-1] = bases[i];
      // printf("right neighbour of %s is %s\n", kmer, right_kmer);
      if(bloom_check(bloom, right_kmer, kmer_size))
      {
        // printf("right neighbour %s of %s kmer present\n", right_kmer, kmer);
        if(contains_left) return 1;
        contains_right = 1;
        // if(bloom_check(edge_bloom, kmer, kmer_size)) return 1;
        // printf("right neighbour %s of %s kmer present\n", tmp_kmer, kmer);
        // break;
      }
    // }
    // if(!contains_left)
    // {
      left_kmer[0] = bases[i];
      // printf("left neighbour of %s is %s\n", kmer, left_kmer);
      if(bloom_check(bloom, left_kmer, kmer_size))
      {
        // printf("left neighbour of %s is %s\n", kmer, left_kmer);
        if(contains_right) return 1;
        contains_left = 1;
        // if(bloom_check(edge_bloom, kmer, kmer_size)) return 1;
        // printf("right neighbour %s of %s kmer present\n", tmp_kmer, kmer);
        // break;
      }
    // }
  }
    // else
    // {
    //   printf("right neighbour %s of %s kmer not present\n", tmp_kmer, kmer);
    // }
  if(contains_left || contains_right)
  {
    // printf("left or right\n");
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
  // uint8_t check;
  uint8_t cnt = 0;
  unsigned int line_cnt = 0, l;

  // printf("\n\n******** Sparse FASTA ********\n\n");

  // memset(kmer, 0, sizeof(kmer));
  kmer[kmer_size] = '\0';
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
    i=0;
    while(i < len - 1)
    {
      if(((strlen(line + i) - 1 < kmer_size)
          || (strlen(line + i) + strlen(last_kmer) - 1 < kmer_size))
          && first)
      {
        // snprintf(last_kmer, strlen(line), "%s", line + i);
        for(l = 0; l <= s; l++)
        {
          // printf("line %s\n", line + i - l);
          if(strlen(line + i - l) - 1 == kmer_size)
          {
            snprintf(last_kmer, kmer_size + 1, "%s", line + i - l);
            break;
          }
        }
        cnt = kmer_size - strlen(line + i) + 1;
        // printf("cnt is %d\n", cnt);
        first = 0;
        break;
      }
      // memset(kmer, 0, sizeof(kmer));
      if(first)
      {
        snprintf(kmer, kmer_size+1, "%s", line+i);
        // if(((s == 1) && !(cnt % 2)) || ((s != 1) && !(cnt % (s + 1))))
        // {
          bloom_add(sparse_bloom, kmer, kmer_size);
          // printf("in first kmer to add: %s\n", kmer);
        // }
        // cnt++;
        // if(cnt == s+1)
        // {
        //   // if(s == 1)
        //   // {
        //     cnt = 0;
        //   // }
        //   // else
        //   // {
        //   //   cnt = 1;
        //   // }
        // }
        snprintf(last_kmer, kmer_size + 1, "%s", line + i);
        // snprintf(edge_kmer, kmer_size + 1, "%s", kmer);
        i=i+(s+1);
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
        memmove(last_kmer, &(last_kmer[cnt]), strlen(&(last_kmer[cnt])));
        // last_kmer[strlen(last_kmer) - 1] = '\0';
        memset(&last_kmer[strlen(last_kmer)-cnt], 0, cnt);
        cnt = 0;
        while(strlen(last_kmer) > 0)
        {
          // memset(kmer, 0, sizeof(kmer));
          j++;
          snprintf(kmer, strlen(last_kmer) + 1, "%s", last_kmer);
          kmer[strlen(last_kmer)] = '\0';
          strncat(kmer, line, kmer_size - strlen(last_kmer));
          /* in case that last line of input file has less characters than
            kmer_size, break at last k-mer. If break statement is missing then
            the last k-mer won't be of kmer_size length */
          if(strlen(kmer) < kmer_size)
          {
            // i = cnt;
            break;
          }
          if(((s == 1) && !(cnt % 2)) || ((s != 1) && !(cnt % (s + 1))))
          {
            // printf("last_kmer %s\n", last_kmer);
            // printf("in else kmer to add: %s\n", kmer);
            bloom_add(sparse_bloom, kmer, kmer_size);
            i = s+1-strlen(last_kmer);
            // printf("strlen = %d\n", i);
          }
          memmove(last_kmer, &(last_kmer[1]), strlen(&(last_kmer[1])));
          last_kmer[strlen(last_kmer) - 1] = '\0';
          cnt++;
          if(cnt == s+1)
          {
            cnt = 0;
          }
          // snprintf(edge_kmer, kmer_size + 1, "%s", kmer);
          // i = cnt;
        }
        first = 1;
      }
      if(left_edge)
      {
        // printf("left edge kmer to add: %s\n", kmer);
        // check = bloom_add(edge_bloom, kmer, kmer_size);
        left_edge = 0;
      }
    }
    first = 0;
  }
  // printf("right edge kmer to add: %s\n", last_kmer);
  bloom_add(sparse_bloom, last_kmer, kmer_size);
  bloom_add(edge_bloom, last_kmer, kmer_size);
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
  for(unsigned int i = 0; i < 4; i++)
  {
    if(*final_contain)
    {
      break;
    }
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
            // neighbour, dist);
    if(dist == 1){
      contains = bloom_check(sparse_bloom, neighbour, kmer_size);
      *final_contain |= contains;
      if(contains)
      {
        break;
      }
      // if (contains)
      // {
      //   if(left)
      //   {
          // printf("left neighbour %s of %s is present in bloom filter\n",
                  // neighbour, query);
        // }
        // else
        // {
          // printf("right neighbour %s of %s is present in bloom filter\n",
                  // neighbour, query);
        // }
      // }
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
  // char right_neighbout[kmer_size + 1], left_neighbour[kmer_size + 1];
  int contains_left = 0, contains_right = 0, contains = 0;
  uint8_t i;
  if(bloom_check(sparse_bloom, query, kmer_size))
  {
    i = 1;
    contains = 0;
    snprintf(&neighbour[s+1], kmer_size - (s+1) + 1, "%s", query);
    // printf("first naighbour %s of %s for %d\n", neighbour, query, i);
    contains_left = strict_contains_neighbours(query, sparse_bloom, (s+1), 1,
          kmer_size, s, neighbour, &contains);
    snprintf(neighbour, kmer_size - (s+1) + 1, "%s", query+s+1);
    // printf("first naighbour %s of %s for %d\n", neighbour, query, i);
    contains = 0;
    contains_right = strict_contains_neighbours(query, sparse_bloom, (s+1), 0,
          kmer_size, s, neighbour, &contains);
    // printf("left %d right %d\n", contains_left, contains_right);
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
    }
    return 1;
  }
  i = 1;
  contains = 0;
  contains_left = contains_right = 0;
  snprintf(&neighbour[i], kmer_size - i + 1, "%s", query);
  // printf("first naighbour %s of %s for %d\n", neighbour, query, i);
  contains_left = strict_contains_neighbours(query, sparse_bloom, i, 1,
        kmer_size, i, neighbour, &contains);
  snprintf(neighbour, kmer_size - i + 1, "%s", query+i);
  // printf("first naighbour %s of %s for %d\n", neighbour, query, i);
  contains = 0;
  contains_right = strict_contains_neighbours(query, sparse_bloom, i, 0,
        kmer_size, i, neighbour, &contains);
  // printf("left %d right %d\n", contains_left, contains_right);
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
  }
  contains_left = contains_right = 0;
  for (i = 0; i < s; i++)
  {
    // printf("********* i = %d ***********\n", i);
    contains = 0;
    snprintf(&neighbour[i+1], kmer_size - i, "%s", query);
    // printf("first naighbour %s of %s for %d\n", neighbour, query, i);
    contains_left = strict_contains_neighbours(query, sparse_bloom, i+1, 1,
          kmer_size, i+1, neighbour, &contains);
    snprintf(neighbour, kmer_size - (s - (i + 1)), "%s", query+(s - (i)));
    // printf("first naighbour %s of %s for %d\n", neighbour, query, i);
    contains = 0;
    contains_right = strict_contains_neighbours(query, sparse_bloom, s - (i), 0,
        kmer_size, s - (i), neighbour, &contains);
    // printf("****************************\n");
    if (contains_left && contains_right)
    {
      return 1;
    }
    if (contains_left || contains_right)
    {
      if(bloom_check(edge_bloom, query, kmer_size))
      {
        return 1;
      }
    }
  }
  return 0;
}
