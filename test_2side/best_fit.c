#include <best_fit.h>

int best_fit_parse ( FILE *fd, struct bloom * sparse_bloom,
      struct bloom * edge_bloom, size_t kmer_size, uint8_t s )
{
  char kmer[kmer_size + 1];
  char *line = NULL, last_kmer[kmer_size + 1], edge_kmer[kmer_size + 1];
  size_t len = 0, kmer_length, line_length;
  ssize_t read;
  uint8_t left_edge = 0, first = 1, i;
  // uint8_t check;
  uint8_t find_index = 0, idx_flag = 0, last_read = 0;
  unsigned int line_cnt = 0, l, file_offset = 0, cnt = 0;
  uint16_t index[s+1], idx = 0, max;
  uint8_t k = sizeof(index)/sizeof(uint16_t);

  // printf("\n\n******** Sparse FASTA ********\n\n");

  // memset(kmer, 0, sizeof(kmer));
  kmer[kmer_size] = '\0';
  for(i = 0; i < k; i++)
  {
    index[i] = 0;
  }
  while(((read = getline(&line, &len, fd)) != -1) || !last_read)
  {
    if(read == -1)
    {
        for(i = 1; i <k; i++)
        {
          if(max < index[i])
          {
            max = index[i];
            idx = i;
          }
        }
        // if(line_cnt) find_index = 1;
        // if(cnt == 1) file_offset += offset;
        // offset = 0;
        // printf("ftell before %ld %d\n", ftell(fd), offset);
        // file_offset -= offset;
        fseek(fd, file_offset, SEEK_SET);
        // printf("ftell %ld\n", ftell(fd));
        find_index = 0;
        left_edge = 1;
        last_read = 1;
        first = 1;
        continue;
    }
    if(strstr(line, ">"))
    {
      // file_offset += strlen(line);
      // printf("idx in strstr is %d\n", find_index);
      if(line_cnt)
      {
        find_index ^= 1;
        bloom_add(edge_bloom, edge_kmer, kmer_size);
        first = 1;
        // printf("idx is %d\n", find_index);
        left_edge = 1;
        // continue;
      }
      else
      {
        // file_offset += strlen(line);
        left_edge = 1;
        line_cnt++;
        // printf("ftell %ld\n", ftell(fd));
        continue;
      }
      if(!find_index)
      {
        for(i = 1; i < k; i++)
        {
          if(max < index[i])
          {
            max = index[i];
            idx = i;
          }
        }
        // if(line_cnt) find_index = 1;
        // if(cnt == 1) file_offset += offset;
        // offset = 0;
        // printf("ftell before %ld\n", ftell(fd));
        // file_offset += offset;
        fseek(fd, file_offset, SEEK_SET);
        // printf("ftell %ld\n", ftell(fd));
        left_edge = 1;
        continue;
      }
      else
      {
        max = index[0];
        idx = 0;
        idx_flag = 1;
        // file_offset += strlen(line);
        file_offset = ftell(fd);
        // printf("ftell find index %ld\n", ftell(fd));
        // if(line_cnt == 2) file_offset += offset;
        for(i = 0; i < k; i++)
        {
          index[i] = 0;
        }
        idx = 0;
        // find_index = 0;
        // file_offset += strlen(line);
        left_edge = 1;
        continue;
      }
    }
    if(find_index) cnt = 0;
    if(idx_flag)
    {
      i = idx;
      idx_flag = 0;
    }
    else
    {
      i = 0;
    }
    while(i < len - 1)
    {
      if(!find_index)
      {
        line_length = strlen(line + i);
        kmer_length = strlen(last_kmer);
        if(((line_length - 1 < kmer_size)
            || (line_length + kmer_length - 1 < kmer_size))
            && first)
        {
          for(l = 0; l <= s; l++)
          {
            // printf("line %s\n", line + i - l);
            if(strlen(line + i - l) - 1 == kmer_size)
            {
              snprintf(last_kmer, kmer_size + 1, "%s", line + i - l);
              break;
            }
          }
          cnt = kmer_size - line_length + 1;
          // printf("cnt is %d\n", cnt);
          first = 0;
          // snprintf(last_kmer, strlen(line), "%s", line + i);
          break;
        }
        // memset(kmer, 0, sizeof(kmer));
        if(first)
        {
          snprintf(kmer, kmer_size+1, "%s", line+i);
          bloom_add(sparse_bloom, kmer, kmer_size);
          // printf("in first kmer to add: %s\n", kmer);

          snprintf(last_kmer, kmer_size + 1, "%s", line + i);
          // snprintf(edge_kmer, kmer_size + 1, "%s", kmer);
          i=i+(s+1);
        }
        else
        {
          unsigned int j = 0;
          // printf("in else, last_kmer = %s\n", last_kmer);
          line_length = strlen(line);
          if(line_length < kmer_size)
          {
            /* If length of line is < kmer_size (it could be the last line of
               sequence read) don't include new line character '\n' */
            *(line + line_length - 1) = '\0';
          }
          memmove(last_kmer, &(last_kmer[cnt]), strlen(&(last_kmer[cnt])));
          // last_kmer[strlen(last_kmer) - 1] = '\0';
          kmer_length = strlen(last_kmer);
          memset(&last_kmer[kmer_length-cnt], 0, cnt);
          kmer_length = strlen(last_kmer);
          cnt = 0;
          while(kmer_length > 0)
          {
            // memset(kmer, 0, sizeof(kmer));
            j++;
            snprintf(kmer, kmer_length + 1, "%s", last_kmer);
            kmer[kmer_length] = '\0';
            strncat(kmer, line, kmer_size - kmer_length);
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
              i = s+1-kmer_length;
              // printf("strlen = %d\n", i);
            }
            memmove(last_kmer, &(last_kmer[1]), kmer_length-1);
            last_kmer[kmer_length - 1] = '\0';
            kmer_length = kmer_length - 1;
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
      else //find_index == 1
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
        if(first)
        {
          snprintf(kmer, kmer_size+1, "%s", line+i);
          index[cnt % (s + 1)] += bloom_check(sparse_bloom, kmer, kmer_size);
          cnt++;
          // printf("in first check kmer: %s %d idx %d\n", kmer, index[(cnt - 1) % (s + 1)],(cnt - 1) % (s + 1));

          // snprintf(last_kmer, kmer_size + 1, "%s", line + i);
          // snprintf(edge_kmer, kmer_size + 1, "%s", kmer);
          i = i + 1;
        }
        else
        {
          unsigned int j = 0;
          // printf("in else, last_kmer = %s\n", last_kmer);
          line_length = strlen(line);
          if(line_length < kmer_size)
          {
            /* If length of line is < kmer_size (it could be the last line of
               sequence read) don't include new line character '\n' */
            *(line + line_length - 1) = '\0';
          }
          memmove(last_kmer, &(last_kmer[cnt]), strlen(&(last_kmer[cnt])));
          last_kmer[kmer_length - 1] = '\0';
          kmer_length -= 1;
          // memset(&last_kmer[strlen(last_kmer)-cnt], 0, cnt);
          // cnt = 0;
          while(kmer_length > 0)
          {
            // memset(kmer, 0, sizeof(kmer));
            j++;
            snprintf(kmer, kmer_length + 1, "%s", last_kmer);
            kmer[kmer_length] = '\0';
            strncat(kmer, line, kmer_size - kmer_length);
            /* in case that last line of input file has less characters than
              kmer_size, break at last k-mer. If break statement is missing then
              the last k-mer won't be of kmer_size length */
            if(strlen(kmer) < kmer_size)
            {
              // i = cnt;
              break;
            }
            index[cnt % (s+1)] += bloom_check(sparse_bloom, kmer, kmer_size);
            cnt++;
            // printf("in else check kmer: %s %d idx %d\n", kmer, index[(cnt-1)%(s+1)],(cnt-1)%(s+1));
            // i = s+1-strlen(last_kmer);
            memmove(last_kmer, &(last_kmer[1]), kmer_length - 1);
            last_kmer[kmer_length - 1] = '\0';
            kmer_length--;
            // snprintf(edge_kmer, kmer_size + 1, "%s", kmer);
            // i = cnt;
          }
          first = 1;
        }
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
