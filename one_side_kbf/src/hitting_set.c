/* finding set of k-mers to store in B, approach 2
 * author: Katarina Prgeša */


#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<bloom.h>
#include<hitting_set.h>
#include<mutate.h>

// returns head of the list of kmers
kmer_node_t *parse_hitting_set(int kmer_size, FILE *f, struct bloom * bloom,
                        struct bloom *edge_bloom)
{

  printf("Parsing fasta format.....\n");

  // create head for set of kmers
  kmer_node_t *head = NULL;
  head = (kmer_node_t *) malloc(sizeof(kmer_node_t));
  if (head == NULL){
    printf("Failed do allocate memory for head node!\n");
    return 0;
  }

  ssize_t read;
  size_t len = 0;
  char *line = NULL;

  char sequence[kmer_size];
  char next_sequence[kmer_size];
  char previous_sequence[kmer_size];

  int first_line = 0;
  int cnt;
  int line_cnt = 0;

  int added = 0;
  // char no_left[kmer_size];
  int overflow = 0;
  int end_flag = 0;

  char last_prev[kmer_size];
  char last[kmer_size];

  memset( sequence, 0, (kmer_size+1)*sizeof(char) );
  memset( next_sequence, 0, (kmer_size+1)*sizeof(char) );
  memset( previous_sequence, 0, (kmer_size+1)*sizeof(char) );

  memset( last_prev, 0, (kmer_size+1)*sizeof(char) );
  memset( last, 0, (kmer_size+1)*sizeof(char) );

  char *no_left = "NNNNNNNNNNNNNNNNNNNN";

  // read line by line
  while( ( read = getline( &line, &len, f ) ) != -1 ){

    line[strcspn( line, "\n" )] = 0; // remove end of line char
    cnt = 0; // counter of read characters

    if ( strchr(line, '>') ){
      first_line = 1;
      if (line_cnt){
        end_flag = 1;
        snprintf(last_prev, kmer_size + 1, "%s", previous_sequence);
        snprintf(last, kmer_size + 1, "%s", sequence);
      }
      continue;
    }
    else{
      line_cnt++;
    }

    // read first kmer from each read
    if ( first_line ){
      first_line = 0;
      // read first kmer - this kmer has no left neighbour
      snprintf( sequence, kmer_size + 1, "%s", line );
      cnt += kmer_size;

      *line++;
      snprintf( next_sequence, kmer_size + 1, "%s", line); // right neighbour
      cnt++;

      line = line + kmer_size;

      // add sequencies to head
      if(line_cnt == 1){

        head->current_kmer = malloc(sizeof(char)*(kmer_size+1));
        // head->previous_kmers = malloc(sizeof(char)*(kmer_size+1));
        // head->next_kmers = malloc(sizeof(char)*(kmer_size+1));
        head->next = NULL;

        snprintf(head->current_kmer, kmer_size + 1, "%s", sequence);
        // snprintf(head->next_kmers, kmer_size + 1, "%s", next_sequence);
        head->next_count = 1;
        head->previous_count = 0;

        bloom_add(edge_bloom, sequence, kmer_size);
        bloom_add(bloom, sequence, kmer_size);


      }

      // if its not first kmer from file to add, update list of kmers but do
      // not add left neighbour to node
      else{

        bloom_add(edge_bloom, sequence, kmer_size);
        bloom_add(bloom, sequence, kmer_size);

        add_to_list( sequence, no_left, next_sequence, kmer_size, head);

      }

    }

    if (overflow){

      overflow = 0;

      if(!end_flag){
        end_flag = 0;
        snprintf( next_sequence, kmer_size +1, "%s", &sequence[1]);
        next_sequence[kmer_size -1] = *line++;
        cnt++;
        // bloom_add(bloom, sequence, kmer_size);
        sequence[kmer_size] = '\0';
        printf("before adding to list: %s, size:\n",sequence, sizeof(sequence) );
        add_to_list( sequence, previous_sequence, next_sequence, kmer_size, head);
      }
      else{
        // bloom_add(edge_bloom, last, kmer_size);
        // bloom_add(bloom, last, kmer_size);


        }
        add_to_list(last, last_prev, no_left, kmer_size, head);

    }

    while(1){

      snprintf( previous_sequence, kmer_size + 1, "%s", next_sequence );
      snprintf( sequence, kmer_size + 1, "%s", &next_sequence[1]);
      sequence[kmer_size - 1] = *line++;
      cnt++;
      // check whether next neighbour continues in the next line
      if(cnt == (read - 1)){
        overflow = 1;
        break;
      }

      snprintf( next_sequence, kmer_size + 1, "%s", &sequence[1]);
      next_sequence[kmer_size -1] = *line++;
      cnt++;

      // bloom_add(bloom, sequence, kmer_size);

      added += 1;

      add_to_list( sequence, previous_sequence, next_sequence, kmer_size, head);


      // check whether last char of right neighbour is the last char in line
      if(cnt == (read - 1)){
        break;
      }

    }

    line = line - read + 1;
  }

  // bloom_add(edge_bloom, sequence, kmer_size);

  free(line);
  return head;
}

//===================================================================================
//===================================================================================
//  adding element to the end of the list
//  or updating existing node
void add_to_list(char *kmer, char *left, char *right, int kmer_size, kmer_node_t *head){

  kmer_node_t *current = head;
  char *no_left = "NNNNNNNNNNNNNNNNNNNN";
  int cnt=0;


  while( (current->next != NULL) && (strcmp(current->current_kmer, kmer)) ){

    current = current->next;
    cnt++;
  }

  printf("Elements in list = %d\n" , cnt);
  if (cnt == 31){
    printf("%s\n", kmer );
    printf("%s\n", current->current_kmer );
  }
  // if node for kmer already exists -> update neighbouring kmers if neccessary
  if(!(strcmp(current->current_kmer, kmer))){
    // before reallocation check whether neighbours exist
    // if they dont exist add them to their previous/next set

    // if(strcmp(left, no_left)){
    //   if (!check_presence(current->previous_kmers, left, kmer_size)){
    //     current->previous_count++;
    //     char *new_p = (char *)realloc( current->previous_kmers,((kmer_size+1)*current->previous_count));
    //
    //     if (new_p == NULL){
    //       printf("Failed to reallocate memory for previous kmers in node!\n");
    //     }
    //     current->previous_kmers = new_p;
    //     free(new_p);
    //     snprintf(current->previous_kmers, kmer_size + 1, "%s", left);
    //
    //   }
    // }



    // if (!check_presence(current->next_kmers, right, kmer_size)){
    //   current->next_count++;
    //   char *new_n = (char *)realloc( current->next_kmers,((kmer_size+1)*current->next_count));
    //
    //   if (new_n == NULL){
    //     printf("Failed to reallocate memory for next kmers in node!\n");
    //   }
    //
    //   current->next_kmers = new_n;
    //   free(new_n);
    //   snprintf(current->next_kmers, kmer_size + 1, "%s", right);
    // }



  }

  // create new node in kmer list and store everything
  else{

    current->next = (kmer_node_t *) malloc(sizeof(kmer_node_t));

    if (  current->next == NULL){
      printf("Failed to allocate memory for new node!\n" );
    }
    current->next->next = NULL;

    current->next->current_kmer = malloc(sizeof(char)*(kmer_size+1));
    // current->next->next_kmers = malloc(sizeof(char)*(kmer_size+1));

    if (  current->next->current_kmer == NULL){
      printf("Failed to allocate memory for current kmer in new node!\n" );
    }

    // if (  current->next->next_kmers == NULL){
    //   printf("Failed to allocate memory for next kmers in new node!\n" );
    // }

    snprintf(current->next->current_kmer, kmer_size + 1, "%s", kmer);
    // snprintf(current->next->next_kmers, kmer_size + 1, "%s", right);

    current->next->next_count = 1;

    // current->next->previous_kmers = malloc(sizeof(char)*(kmer_size+1));

    // if (  current->next->previous_kmers == NULL){
    //   printf("Failed to allocate memory for previous kmers in new node!\n" );
    // }

    // if(strcmp(left, no_left)){
    //   snprintf(current->next->previous_kmers, kmer_size + 1, "%s", left);
    //   current->next->previous_count = 1;
    // }

    // add kmer to the set of previous kmers
    // char *new = (char *)realloc( current->next->previous_kmers,((kmer_size+1)*current->next->previous_count));

    // if (  new == NULL){
    //   printf("Failed to reallocate memory for previous kmers in node!\n" );
    // }

    // current->next->previous_kmers = new;
    // snprintf(current->next->previous_kmers, kmer_size + 1, "%s", kmer);

    // current->next->previous_count++;
  }

  return;

}

//==============================================================================
//==============================================================================
// check whether data_set contains query of length query_length
// if contains return 1, else return 0
int check_presence(char *data_set, char *query, int query_length){

  char tmp[query_length];
  memset(tmp, 0, query_length*sizeof(char));

  for (int i = 0; data_set[i] != NULL; i += query_length ){
    snprintf(tmp, query_length + 1, "%s", &data_set[i]);
    if (!strcmp(tmp, query)){
      return 1;
    }
  }

  return 0;
}
