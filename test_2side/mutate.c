#include <mutate.h>

int mutate ( FILE *fd, size_t kmer_size, uint8_t mutate )
{
  FILE *new_fd;
  uint8_t rnd_idx;
  unsigned int queries = 0, mutated = 0;
  char *line = NULL, kmer[kmer_size + 2];
  size_t len = 0;
  char bases[] = {'A', 'C', 'G', 'T'}, base, *base_idx;
  unsigned int counter = 0, cnts;
  srand(time(0));
  new_fd = fopen("../test_files/mutate.txt", "w");
  memset(kmer, '1', sizeof(kmer));
  cnts = rand() % (5) + 1;
  base = bases[rand() % sizeof(bases)];
  while(queries < QUERIES)
  {
    if(getline(&line, &len, fd) == -1)
    {
      memset(kmer, '1', sizeof(kmer));
      fseek(fd, 0, SEEK_SET);
      getline(&line, &len, fd);
    }
    if(strstr(line, ">"))
    {
      memset(kmer, '1', sizeof(kmer));
      getline(&line, &len, fd);
    }
    if(strlen(line) - 1 < kmer_size)
    {
      memset(kmer, '1', sizeof(kmer));
      fseek(fd, 0, SEEK_SET);
      continue;
    }
    if(strlen(kmer) < kmer_size + 1)
    {
      // kmer[strlen(kmer)] = '\0';
      strncat(kmer, line, sizeof(kmer) - strlen(kmer) - 2);
      kmer[strlen(kmer)] = '\n';
      if( (counter >= cnts) && mutate )
      {
        cnts = rand() % (5) + 1;
        counter = 0;
        base_idx = strchr(kmer, base);
        while(base_idx)
        {
          *base_idx = bases[rand() % sizeof(bases)];
          base_idx = strchr(base_idx, base);
        }
        mutated++;
        base = bases[rand() % sizeof(bases)];
      }
      fputs(kmer, new_fd);
      // printf("kmer in if is %s",  kmer);
      queries++;
      if(queries == QUERIES)
      {
        break;
      }
      fseek(new_fd, 0, SEEK_END);
    }
    memset(kmer, 0, sizeof(kmer));
    rnd_idx = rand() % (strlen(line) - 1);
    if(strlen(line + rnd_idx) - 1 < kmer_size)
    {
      // strncat(kmer, line + rnd_idx, strlen(line + rnd_idx) - 1);
      snprintf(kmer, strlen(line + rnd_idx), "%s", line + rnd_idx);
    }
    else
    {
      // strncat(kmer, line + rnd_idx, kmer_size);
      snprintf(kmer, kmer_size + 1, "%s", line + rnd_idx);
      kmer[strlen(kmer)] = '\n';
      if( (counter >= cnts) && mutate )
      {
        cnts = rand() % (5) + 1;
        counter = 0;
        base_idx = strchr(kmer, base);
        while(base_idx)
        {
          *base_idx = bases[rand() % sizeof(bases)];
          base_idx = strchr(base_idx, base);
        }
        mutated++;
        base = bases[rand() % sizeof(bases)];
      }
      fputs(kmer, new_fd);
      // printf("kmer is %s",  kmer);
      queries++;
      fseek(new_fd, 0, SEEK_END);
    }
    counter++;
  }
  if(line)
  {
    free(line);
  }
  fclose(new_fd);
  return mutated;
}
