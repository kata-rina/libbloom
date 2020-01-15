#include<hitting_set.h>
#include<relaxed_sparse.h>
#include<bloom.h>
#include<stdio.h>
#include<stdlib.h>
#include<stdint.h>

//==============================================================================
//==============================================================================

int decide_present(char *query, int size, int contains_l, int contains_r,
                    struct bloom *edge_bloom){

  if ( contains_l && contains_r)
    return 1;

  if ( contains_l || contains_r){
    if (bloom_check(edge_bloom, query, size) == 1){
      return 1;
    }
  }
  return 0;
}


// //==============================================================================
// //==============================================================================

int relaxed_contains(char *query, int size, int s,
                      struct bloom *bloom, struct bloom *edge_bloom ){

    if (bloom_check(bloom, query, size) == 1){

      if (relaxed_contains_neighbours(query, (s + 1), (s + 1), bloom, edge_bloom, size)){
        return 1;
      }

    }

    if(relaxed_contains_neighbours(query,  1,  1, bloom, edge_bloom, size))
      return 1;

    for(int i = 0; i < s + 1; i++){

      if (relaxed_contains_neighbours(query, i, (s - (i)), bloom, edge_bloom, size)){
        return 1;
      }

    }



    return 0;
}



//==============================================================================
//==============================================================================

int relaxed_contains_neighbours(char *query, int left_dist, int right_dist,
            struct bloom * bloom, struct bloom *edge_bloom, int kmer_size)
{
  int contains_left = 0;
  int contains_right = 0;

  left_neighbours_relaxed(query, bloom, 1, kmer_size, left_dist, &contains_left);
  right_neighbours_relaxed(query, bloom, 1, kmer_size, right_dist, &contains_right);

  return decide_present(query, kmer_size, contains_left, contains_right, edge_bloom);
}

//==============================================================================
//==============================================================================

void left_neighbours_relaxed ( char * query,
      struct bloom * sparse_bloom, int dist, int kmer_size,
      int s, int *contains)
{

  int contains_left = 0;
  if ( s == 0){
    return;
  }
  if (dist > s){
    return;
  }


  else{
    char bases[] = {'A', 'T', 'G', 'C'};
    char neighbour[kmer_size+1];
    snprintf(&neighbour[1], kmer_size, "%s", query);

    for (int i = 0; i < 4; i++){
      neighbour[0] = bases[i];
      neighbour[kmer_size] = NULL;
      contains_left = bloom_check(sparse_bloom, neighbour, kmer_size);
      if (contains_left || *contains){
        *contains |= contains_left;
        return;
      }
      left_neighbours_relaxed(neighbour, sparse_bloom, (dist + 1), kmer_size, s, contains);
    }
  }
  *contains |= contains_left;
  return;

  }


//==============================================================================
//==============================================================================

void right_neighbours_relaxed ( char * query,
      struct bloom * sparse_bloom, int dist, int kmer_size,
      int s, int *contains)
{

  int contains_right = 0;
  if (s == 0){
    return;
  }
  if (dist > s){
    return;
  }


  else{
    char bases[] = {'A', 'T', 'G', 'C'};
    char neighbour[kmer_size+1];

    snprintf(neighbour, kmer_size, "%s", &query[1]);

    for (int i = 0; i < 4; i++){
      neighbour[kmer_size - 1] = bases[i];
      neighbour[kmer_size] = NULL;
      contains_right = bloom_check(sparse_bloom, neighbour, kmer_size);

      if(contains_right || *contains){
        *contains |= contains_right;
        return;
      }

      right_neighbours_relaxed(neighbour, sparse_bloom, (dist + 1), kmer_size, s, contains);
    }

  }

  *contains |= contains_right;
  return;

}
