/* finding set of k-mers to store in B, approach 2
 * author: Katarina Prgeša */


#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<bloom.h>
#include<hitting_set.h>

// returns head of the list of kmers
kmer_node_t *parse_hitting_set(int kmer_size, int skip_length, FILE *f,struct bloom * bloom){

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
  char sequence[kmer_size];
  char next_sequence[kmer_size];
  char previous_sequence[kmer_size];
  int first_line = 0;
  int cnt;
  int line_cnt = 0;
  char *line = NULL;

  char no_left[kmer_size];
  int overflow = 0;
  int n = 0;
  int end_flag = 0;

  memset( sequence, 0, (kmer_size+1)*sizeof(char) );
  memset( next_sequence, 0, (kmer_size+1)*sizeof(char) );
  memset( previous_sequence, 0, (kmer_size+1)*sizeof(char));

  // read line by line
  while( ( read = getline( &line, &len, f ) ) != -1 ){

    line[strcspn( line, "\n" )] = 0; // remove end of line char
    cnt = 0; // counter of read characters

    if ( strchr(line, '>') ){
      first_line = 1;
      if (line_cnt){
        end_flag = 1;
        // tu treba spremiti sequence i previous sequence da se ne prepišu
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
        head->previous_kmers = malloc(sizeof(char)*(kmer_size+1));
        head->next_kmers = malloc(sizeof(char)*(kmer_size+1));
        head->next = NULL;

        snprintf(head->current_kmer, kmer_size + 1, "%s", sequence);
        snprintf(head->next_kmers, kmer_size + 1, "%s", next_sequence);
        head->next_count = 1;
        head->previous_count = 0;

      }

      // if its not first kmer from file to add, update list of kmers but do
      // not add left neighbour to node
      else{

        add_to_list( sequence, no_left, next_sequence, kmer_size, head);

      }

    }

    if (overflow){

      printf("%d overflow occured\n", line_cnt );
      overflow = 0;

      if(!end_flag){
        end_flag = 0;
        snprintf( next_sequence, kmer_size +1, "%s", &sequence[1]);
        next_sequence[kmer_size -1] = *line++;
        cnt++;
        n++;
        printf("%d)%s size %d\n", n, previous_sequence, sizeof(previous_sequence));
        printf("%d)%s size %d\n", n, sequence, sizeof(sequence));
        printf("%d)%s size %d\n", n, next_sequence, sizeof(next_sequence));
        add_to_list( sequence, previous_sequence, next_sequence, kmer_size, head);
      }
      else{
      //krive sekvence spremaš!!
        printf("adding to list with no right neighbour");
        printf("%d)%s size %d\n", n, previous_sequence, sizeof(previous_sequence));
        printf("%d)%s size %d\n", n, sequence, sizeof(sequence));
        add_to_list(sequence, previous_sequence, no_left, kmer_size, head);
      }
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

      n++;
      // printf("%d)%s\n", n, previous_sequence);
      // printf("%d)%s\n", n, sequence);
      // printf("%d)%s\n", n, next_sequence);
      add_to_list( sequence, previous_sequence, next_sequence, kmer_size, head);


      // check whether last char of right neighbour is the last char in line
      if(cnt == (read - 1)){
        break;
      }

    }

    line = line - read + 1;
  }
  // free(line);
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
  // if node for kmer already exists -> update neighbouring kmers if neccessary
  if(!(strcmp(current->current_kmer, kmer))){
    printf("već postojim\n" );
    // before reallocation check whether neighbours exist
    // if they dont exist add them to their previous/next set
    if(strcmp(left, no_left)){
      if (!check_presence(current->previous_kmers, left, kmer_size)){
        current->previous_count++;
        char *new_p = (char *)realloc( current->previous_kmers,((kmer_size+1)*current->previous_count));

        if (new_p == NULL){
          printf("realokacija za lijeve susjede nije uspjela\n");
        }
        current->previous_kmers = new_p;
        free(new_p);
        snprintf(current->previous_kmers, kmer_size + 1, "%s", left);

      }
    }

    if (!check_presence(current->next_kmers, right, kmer_size)){
      current->next_count++;
      char *new_n = (char *)realloc( current->next_kmers,((kmer_size+1)*current->next_count));

      if (new_n == NULL){
        printf("realokacija za desne susjede nije uspjela\n");
      }

      current->next_kmers = new_n;
      free(new_n);
      snprintf(current->next_kmers, kmer_size + 1, "%s", right);
    }

  }

  // create new node in kmer list and store everything
  else{

    current->next = (kmer_node_t *) malloc(sizeof(kmer_node_t));

    if (  current->next == NULL){
      printf("alokacija za novi element liste nije uspjela\n" );
    }
    current->next->next = NULL;

    current->next->current_kmer = malloc(sizeof(char)*(kmer_size+1));
    current->next->next_kmers = malloc(sizeof(char)*(kmer_size+1));

    if (  current->next->current_kmer == NULL){
      printf("alokacija za current kmer liste nije uspjela\n" );
    }

    if (  current->next->next_kmers == NULL){
      printf("alokacija za next kmer liste nije uspjela\n" );
    }

    snprintf(current->next->current_kmer, kmer_size + 1, "%s", kmer);
    snprintf(current->next->next_kmers, kmer_size + 1, "%s", right);

    current->next->next_count = 1;

    current->next->previous_kmers = malloc(sizeof(char)*(kmer_size+1));

    if (  current->next->previous_kmers == NULL){
      printf("alokacija za previous kmer liste nije uspjela\n" );
    }

    if(strcmp(left, no_left)){
      snprintf(current->next->previous_kmers, kmer_size + 1, "%s", left);
      current->next->previous_count = 1;
    }

    // add kmer to the set of previous kmers
    char *new = (char *)realloc( current->next->previous_kmers,((kmer_size+1)*current->next->previous_count));

    if (  new == NULL){
      printf("2realokacija za previous kmer liste nije uspjela\n" );
    }

    current->next->previous_kmers = new;
    snprintf(current->next->previous_kmers, kmer_size + 1, "%s", kmer);

    current->next->previous_count++;
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
