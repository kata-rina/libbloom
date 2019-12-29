/*********** One sided k-mer bloom filter implementation***********************/
/* author: Katarina Prgeša
 *
 *
 *
 *
 *
 *
 */
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<bloom.h>
#include<onesidekbf.h>

/* Function adds k-mers in bloomfilter */
int genom_add(FILE * fp, int kmer_size, struct bloom * bloom){

  ssize_t read;
  size_t len = 0;
  char *line = NULL;
  char sequence[kmer_size];
  char previous_sequence[kmer_size];

  memset( sequence, 0, kmer_size*sizeof(char) );
  memset( previous_sequence, 0, kmer_size*sizeof(char) );

  int first_line = 1;
  int cnt;

  // read line by line
  while( ( read = getline( &line, &len, fp ) ) != -1 ){

    line[strcspn( line, "\n" )] = 0; // remove end of line char
    cnt = 0; // counter of read characters

    // read first kmer from file
    if ( first_line ){
      first_line = 0;

      snprintf( sequence, kmer_size + 1, "%s", line );

      for ( int i=0; i<kmer_size; i++ ){
        *line++;
      }

      bloom_add( bloom, sequence, kmer_size );

      cnt += kmer_size;

      snprintf( previous_sequence, kmer_size + 1, "%s", sequence );

    }

    while (cnt < read - 1){

      for ( int i = 0; i < kmer_size - 1; i++ ){
        sequence[i] = previous_sequence[i+1];
      }

      sequence[ kmer_size - 1 ] = *line++;

      bloom_add( bloom, sequence, kmer_size );
      // printf("%s\n", sequence);
      cnt++;

      for ( int i = 0; i < kmer_size; i++ ){
        previous_sequence[i] = sequence[i];
      }

    }

    // restore pointer to line
    line = line - read + 1;
  }

  free( line );
  return 1;
  }

/* This function checks whether bloom filter contains kmer or not */
int onesided_kbf_check( char *kmer, int kmer_size, struct bloom * bloom ){

      if ( bloom_check( bloom, kmer, kmer_size) == 1 ){

        int i;
        char dna_base[] = { 'A', 'T', 'G', 'C'};

        char left_neighbour[kmer_size];
        memset( left_neighbour, 0, kmer_size*sizeof(char) );
        snprintf( left_neighbour, kmer_size + 1, "%s", kmer);

        for (i = kmer_size - 1; i > 0; i--){
          left_neighbour[i] = left_neighbour[i-1];
        }

        char right_neighbour[kmer_size];
        memset( right_neighbour, 0, kmer_size*sizeof(char) );
        snprintf( right_neighbour, kmer_size + 1, "%s", kmer);

        for (i = 0; i < kmer_size; i++){
          right_neighbour[i] = right_neighbour[i+1];
        }


        // check for presence of neighbouring kmers
        for ( i = 0; i < 4; i++ ){
            left_neighbour[0] = dna_base[i];
            right_neighbour[kmer_size - 1] = dna_base[i];

            // check for presence of the left neighbour
            if ( bloom_check( bloom, left_neighbour, kmer_size) == 1 ) {return 1;}
            if ( bloom_check( bloom, right_neighbour, kmer_size) == 1) {return 1;}
          }
        }

      return 0;
}
