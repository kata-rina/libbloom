#include "bloom.h"
#include <stdio.h>
#include <string.h>

int main(void){

  struct bloom bloom_filter;
  char *buffer = "Katarina";
  int buff_len;

  buff_len = strlen(buffer);

  bloom_init(&bloom_filter, 10000, 0.01);
  bloom_add(&bloom_filter, buffer, buff_len);

  if (bloom_check(&bloom_filter, buffer, buff_len)){
    printf( "I may be there :D\n" );
  }

  char *buff = "Kata";
  buff_len = strlen(buff);

  if (bloom_check(&bloom_filter, buff, buff_len)){
    printf( "I may be there :D\n" );
  }
  else{
    printf("Definately not in filter!\n");
  }

  return 0;

}
