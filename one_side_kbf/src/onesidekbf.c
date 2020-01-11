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
#include<mutate.h>

/* Function parses fasta file and adds k-mers in bloomfilter */
int parse_fasta(FILE * fp, int kmer_size, struct bloom * bloom){

  ssize_t read;
  size_t len = 0;
  char *line = NULL;
  char sequence[kmer_size];
  char previous_sequence[kmer_size];
  int first_line;
  int cnt;

  int added=0; // counter of added kmers in bloom filter

  printf("Parsing fasta format and storing to bloom filter......\n");

  memset( sequence, 0, kmer_size*sizeof(char) );
  memset( previous_sequence, 0, kmer_size*sizeof(char) );



  // read line by line
  while( ( read = getline( &line, &len, fp ) ) != -1 ){

    line[strcspn( line, "\n" )] = 0; // remove end of line char
    cnt = 0; // counter of read characters

    if ( strchr(line, '>') ){
      first_line = 1;
      continue;}

    // read first kmer from each read
    if ( first_line ){
      first_line = 0;

      snprintf( sequence, kmer_size + 1, "%s", line );

      for ( int i=0; i<kmer_size; i++ ){
        *line++;
      }

      if (bloom_add( bloom, sequence, kmer_size ) == 0) { added++;}

      cnt += kmer_size;

      snprintf( previous_sequence, kmer_size + 1, "%s", sequence );

    }

    while (cnt < read - 1){

      for ( int i = 0; i < kmer_size - 1; i++ ){
        sequence[i] = previous_sequence[i+1];
      }

      sequence[ kmer_size - 1 ] = *line++;

      if (bloom_add( bloom, sequence, kmer_size ) == 0) {added++;}
      cnt++;

      for ( int i = 0; i < kmer_size; i++ ){
        previous_sequence[i] = sequence[i];
      }

    }

    // restore pointer to line
    line = line - read + 1;
    // if (added == QUERIES) {
    //
    //   break;}
  }

  printf("* %d kmers added to bloom filter\n\n", added);
  free( line );
  return 1;
  }

/* This function checks whether one sided bloom filter contains kmer or not */
int onesided_kbf_contains( char *kmer, int kmer_size, struct bloom * bloom ){

      if ( bloom_check( bloom, kmer, kmer_size) == 1 ){

        int i;
        char dna_base[] = { 'A', 'T', 'G', 'C'};

        char left_neighbour[kmer_size];
        memset( left_neighbour, 0, kmer_size*sizeof(char) );
        snprintf( &left_neighbour[1], kmer_size, "%s", kmer);


        char right_neighbour[kmer_size];
        memset( right_neighbour, 0, kmer_size*sizeof(char) );
        snprintf( right_neighbour, kmer_size, "%s", &kmer[1]);


        // check for presence of neighbouring kmers
        for ( i = 0; i < 4; i++ ){
            left_neighbour[0] = dna_base[i];
            right_neighbour[kmer_size - 1] = dna_base[i];

            if ( bloom_check( bloom, left_neighbour, kmer_size) == 1 ) {return 1;}
            if ( bloom_check( bloom, right_neighbour, kmer_size) == 1) {return 1;}
          }
        }

      return 0;
}
