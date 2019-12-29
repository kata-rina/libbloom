#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<bloom.h>
#include<onesidekbf.h>

int main(void){

  // file to read from
  FILE *f;
  char filename[] = "/mnt/Jupiter/FAKS/Diplomski/3_semestar/Bioinformarika/genom/test.txt";
  f = fopen(filename, "r");

  // bloom filter to store kmers
  struct bloom bloom_filter;
  bloom_init(&bloom_filter, 10000000, 0.01);

  char *buffer = "TCAACGGGCGGATATCTTGATCAGCTTGGCATAGTAGAAGCCAT";

  // add kmers to bloom filter
  genom_add(f, strlen(buffer), &bloom_filter);
  printf(" Cijeli genom je obraden\n");
  fclose(f);


  if (onesided_kbf_check(buffer,  strlen(buffer), &bloom_filter)){
    printf( "I may be there :D\n" );
  }
  return 0;
}
