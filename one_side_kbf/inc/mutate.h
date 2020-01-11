#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>

#define QUERIES           5600000

int mutate ( FILE *fd, size_t kmer_size, uint8_t mutate );
