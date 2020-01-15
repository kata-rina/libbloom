/* author: Katarina Prgeša */

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<bloom.h>
#include<onesidekbf.h>
#include<mutate.h>
#include<sys/time.h>
#include<hitting_set.h>
#include<relaxed_sparse.h>
#include<stdint.h>


int test_kbf(struct bloom * bloom_filter, FILE *f, int kmer_size, int type);
float calculate_fpr(int queries, int positive, int mutated);
int test_relaxed(struct bloom *bloom, struct bloom *edge, int kmer_size, FILE *f, int s);



int main(void){

  //==============================================================================
  //==============================================================================
  //==============================================================================

  // // for time stats
  struct timeval tstart, tstop;
  //
  // file to read from
  FILE *f;
  int fsize;
  // char filename[] = "/mnt/Jupiter/FAKS/Diplomski/3_semestar/Bioinformarika/genom/GCF_000006765.1_ASM676v1_genomic.fna";
  char filename[] = "/mnt/Jupiter/FAKS/Diplomski/3_semestar/Bioinformarika/genom/sparse.txt";
  f = fopen(filename, "r");
  //
  fseek(f, 0, SEEK_END);
  fsize = ftell(f);
  fseek(f, 0, SEEK_SET);

  if (fsize < 10000){
    fsize *=1000;
  }
  else{
    fsize *= 4;
  }
  //
  // init bloom filter to store kmers
  struct bloom bloom_filter;
  bloom_init(&bloom_filter, fsize, 0.29);
  //
  // parse fasta and add kmers to bloom filter
  gettimeofday(&tstart, NULL);
  parse_fasta(f, KMER_SIZE, &bloom_filter);
  gettimeofday(&tstop, NULL);
  printf("* Time needed to store kmers to bloom filter = %f s\n\n",
            (double) (tstop.tv_usec - tstart.tv_usec) / 1000000 +
            (double) (tstop.tv_sec - tstart.tv_sec));
  //
  //
  // mutate k-mers from original file for testing
  int mutated_nmr = mutate(f, KMER_SIZE, 1);
  int non_mutated_nmr = mutate(f, KMER_SIZE, 0);



  printf("\n" );
  printf("********** One sided kbf performance testing with mutation **********\n" );

  // store with mutation
  FILE *f_mutated;
  f_mutated = fopen("/mnt/Jupiter/FAKS/Diplomski/3_semestar/Bioinformarika/genom/mutate.txt", "r");

  gettimeofday(&tstart, NULL);
  int match_mut = test_kbf(&bloom_filter, f_mutated , KMER_SIZE, 1);
  gettimeofday(&tstop, NULL);

  printf("* Number of mutated k-mers => %d\n", mutated_nmr );
  printf("* Positive queries => %d \n", match_mut );
  printf("* Operating time = %f s\n ",
            (double) (tstop.tv_usec - tstart.tv_usec) / 1000000 +
            (double) (tstop.tv_sec - tstart.tv_sec));

  float fpr = calculate_fpr(QUERIES, match_mut, mutated_nmr);
  printf("* False positive rate => %f%s\n", fpr, "%");




  printf("\n" );
  printf("********** One sided kbf performance testing without mutation **********\n" );

  //store without mutation
  FILE *f_nonmutated;
  f_nonmutated = fopen("/mnt/Jupiter/FAKS/Diplomski/3_semestar/Bioinformarika/genom/nonmutate.txt", "r");

  gettimeofday(&tstart, NULL);
  int match = test_kbf(&bloom_filter, f_nonmutated , KMER_SIZE, 1);
  gettimeofday(&tstop, NULL);

  printf("* Number of mutated k-mers => %d\n", non_mutated_nmr );
  printf("* Positive queries => %d \n", match );
  printf("* Operating time = %f s\n ",
            (double) (tstop.tv_usec - tstart.tv_usec) / 1000000 +
            (double) (tstop.tv_sec - tstart.tv_sec));

  //restoring pointers to begginig of files
  fseek(f_mutated, 0, SEEK_SET);
  fseek(f_nonmutated, 0, SEEK_SET);




  printf("\n" );
  printf("********** Regular kbf performance testing with mutation **********\n" );

  gettimeofday(&tstart, NULL);
  match_mut = test_kbf(&bloom_filter, f_mutated , KMER_SIZE, 0);
  gettimeofday(&tstop, NULL);

  printf("* Number of mutated k-mers => %d\n", mutated_nmr );
  printf("* Positive queries => %d \n", match_mut );
  printf("* Operating time = %f s\n ",
            (double) (tstop.tv_usec - tstart.tv_usec) / 1000000 +
            (double) (tstop.tv_sec - tstart.tv_sec));

  fpr = calculate_fpr(QUERIES, match_mut, mutated_nmr);
  printf("* False positive rate => %f%s\n", fpr, "%");
  //
  //
  printf("\n" );
  printf("********** Regular kbf performance testing without mutation **********\n" );

  gettimeofday(&tstart, NULL);
  match = test_kbf(&bloom_filter, f_nonmutated , KMER_SIZE, 0);
  gettimeofday(&tstop, NULL);

  printf("* Number of mutated k-mers => %d\n", non_mutated_nmr );
  printf("* Positive queries => %d \n", match );
  printf("* Operating time = %f s\n ",
            (double) (tstop.tv_usec - tstart.tv_usec) / 1000000 +
            (double) (tstop.tv_sec - tstart.tv_sec));


  // fclose(f_nonmutated);
  // fclose(f_mutated);
  // fclose(f);

  //=================================================================================
  //=================================================================================
  //=================================================================================

  printf("\n" );
  printf("********** Relaxed kbf performance testing mutation **********\n" );


  struct bloom bloom;
  bloom_init(&bloom, fsize, 0.255);

  struct bloom edge_bloom;
  bloom_init(&edge_bloom, 20000, 0.255);

  fseek(f, 0, SEEK_SET);

  gettimeofday(&tstart, NULL);
  parse_hitting_set(KMER_SIZE, f, &bloom, &edge_bloom);
  gettimeofday(&tstop, NULL);


  printf("* Operating time = %f s\n ",
            (double) (tstop.tv_usec - tstart.tv_usec) / 1000000 +
            (double) (tstop.tv_sec - tstart.tv_sec));

  fseek(f_mutated, 0, SEEK_SET);
  fseek(f_nonmutated, 0, SEEK_SET);
  int s = 1;
  int kmer_size = 20;

  int a = test_relaxed( &bloom, &edge_bloom, kmer_size, f_mutated, s );

  fpr = calculate_fpr(QUERIES, a, mutated_nmr);

  printf("False positive rate => %f\n", fpr );
  //
  // gettimeofday(&tstop, NULL);
  //
  printf("%d\n", a);
  //
  // fclose(f);

  return 0;
}




// function for testing relaxed_sparse kbf performance
int test_relaxed(struct bloom *bloom, struct bloom *edge, int kmer_size, FILE *f, int s){


  ssize_t read;
  size_t len = 0;
  char *line = NULL;
  char sequence[kmer_size];
  int found = 0; // counter of positive queries
  int cnt = 0;

  // bloom_print(edge);
  // bloom_print(bloom);



  while( ( read = getline( &line, &len, f ) ) != -1 ){

    *(line + kmer_size) = '\0';

    if (relaxed_contains( line, kmer_size, s, bloom, edge ))
    {
      found++;
    }
    else{
      // printf("%s\n", line );
    }
  }

  free( line );
  return found;

}



// function for testing kmer bloom filter performance
// if type => 1 - test one sided KBF
// if type => 0 - test regular KBF
int test_kbf(struct bloom * bloom_filter, FILE *f, int kmer_size, int type){

  ssize_t read;
  size_t len = 0;
  char *line = NULL;
  char sequence[kmer_size];
  int found = 0; // counter of positive queries

  memset( sequence, 0, kmer_size*sizeof(char) );

  while( ( read = getline( &line, &len, f ) ) != -1 ){

    snprintf( sequence, kmer_size + 1, "%s", line );

    if(type == 1){
      if (onesided_kbf_contains( sequence, kmer_size, bloom_filter )){found++;}
    }
    else{
      if (bloom_check(bloom_filter, sequence, kmer_size ) == 1){found++;}
    }

  }

  free( line );
  return found;

}

// returns false positive rate in percantage
float calculate_fpr(int queries, int positive, int mutated){

  float fpr;
  fpr = (float)(((positive + mutated) - queries) / (float)queries);
  fpr *= 100;
  return fpr;
}
