#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdbool.h>
#include <math.h>

const double A = 0.618033988;

// ============ DEL 1: HASHTABELL MED TEKSTNØKLER ============

typedef struct Node {
  char *name;
  struct Node *next;
} Node;

typedef struct {
  Node **table;
  int size;
  int count;
  unsigned long long collisions;
} HashTable;

// Hash-funksjon for strenger
unsigned int hash_string(const char *str) {
  unsigned int hash = 0;
  int i = 0;
  while (str[i] != '\0') {
    hash = hash * 31 + (unsigned char)str[i];
    i++;
  }
  return hash;
}

HashTable* create_hash_table(int size) {
  HashTable *ht = malloc(sizeof(HashTable));
  ht->size = size;
  ht->count = 0;
  ht->collisions = 0;
  ht->table = calloc(size, sizeof(Node*));
  return ht;
}

void insert_name(HashTable *ht, const char *name) {
  const unsigned int index = hash_string(name) % ht->size;
  Node *current = ht->table[index];
  while (current != NULL) {
    if (strcmp(current->name, name) == 0) {
      return;
    }
    current = current->next;
  }

  // Create new node
  Node *new_node = malloc(sizeof(Node));
  new_node->name = malloc(strlen(name) + 1);
  strcpy(new_node->name, name);
  new_node->next = ht->table[index];

  // Check for collision
  if (ht->table[index] != NULL) {
    printf("Kollisjon: '%s' og '%s'\n", name, ht->table[index]->name);
    ht->collisions++;
  }

  ht->table[index] = new_node;
  ht->count++;
}

bool search_name(HashTable *ht, const char *name) {
  unsigned int index = hash_string(name) % ht->size;
  Node *current = ht->table[index];
  int probes = 0;

  while (current != NULL) {
    if (probes > 0) {
      printf("Kollisjon ved søk: '%s' og '%s'\n", name, current->name);
    }
    if (strcmp(current->name, name) == 0) {
      return true;
    }
    current = current->next;
    probes++;
  }

  return false;
}

void free_hash_table(HashTable *ht) {
  for (int i = 0; i < ht->size; i++) {
    Node *current = ht->table[i];
    while (current != NULL) {
      Node *temp = current;
      current = current->next;
      free(temp->name);
      free(temp);
    }
  }
  free(ht->table);
  free(ht);
}

void format_tid(char *buf, size_t size, int sec, long nsec) {
    if (sec >= 1.0) {
        snprintf(buf, size, "%.3f s", sec + nsec / 1000000000.0);
    } else if (nsec > 50000000) {
        snprintf(buf, size, "%.1f ms", nsec / 1000000.0);
    } else if (nsec > 50000) {
        snprintf(buf, size, "%.1f µs", nsec / 1000.0);
    } else {
        snprintf(buf, size, "%ld ns", nsec);
    }
}

void del1() {
  printf("\n========== DEL 1: HASHTABELL MED TEKSTNØKLER ==========\n\n");

  FILE *file = fopen("navn.txt", "r");
  if (file == NULL) {
    printf("Kunne ikke åpne navn.txt\n");
    return;
  }

  int name_count = 0;
  char buffer[256];
  while (fgets(buffer, sizeof(buffer), file)) {
    name_count++;
  }
  rewind(file);

  int table_size = (int)(name_count / 0.75);
  if (table_size % 2 == 0) table_size++;

  printf("Antall navn: %d\n", name_count);
  printf("Tabellstørrelse: %d\n\n", table_size);

  HashTable *ht = create_hash_table(table_size);

  while (fgets(buffer, sizeof(buffer), file)) {
    buffer[strcspn(buffer, "\n")] = 0;
    if (strlen(buffer) > 0) {
      insert_name(ht, buffer);
    }
  }
  fclose(file);

  printf("\n--- Statistikk ---\n");
  printf("Totalt antall personer: %d\n", ht->count);
  printf("Totalt antall kollisjoner: %llu\n", ht->collisions);
  printf("Lastfaktor: %.3f\n", (double)ht->count / ht->size);
  printf("Kollisjoner per person: %.3f\n", (double)ht->collisions / ht->count);

  // Test search with a few names
  printf("\n--- Oppslag ---\n");
  const char *test_names[] = {"Ole", "Kari", "Per", "InvalidName", "Yazan Samer,Zarka"};
  for (int i = 0; i < 5; i++) {
    printf("Søker etter '%s': %s\n", test_names[i],
           search_name(ht, test_names[i]) ? "FUNNET" : "IKKE FUNNET");
  }
  free_hash_table(ht);
}

// ============ DEL 2: HASHTABELLER OG YTELSE ============

typedef struct {
  int *table;
  unsigned int size;
  unsigned int count;
  unsigned long long collisions;
} HashTableV2;

unsigned int log2_pow2(const unsigned int size) {
  return (unsigned)log2(size);
}

HashTableV2* create_open_hash_table(const int size) {
  HashTableV2 *ht = malloc(sizeof(HashTableV2));
  ht->size = size;
  ht->count = 0;
  ht->collisions = 0;
  ht->table = malloc(size * sizeof(int));
  for (int i = 0; i < size; i++) {
    ht->table[i] = -1;
  }
  return ht;
}

unsigned int hash1(const int key, const unsigned int m) {
  return key % m;
}

unsigned int hash2(const int key, const unsigned int m) {
  return key % (m - 1) + 1;
}

int insert_linear(HashTableV2 *ht, int key) {
  unsigned int pos = hash1(key, ht->size);
  const unsigned start_pos = pos;

  if (ht->table[pos] == -1) {
    ht->table[pos] = key;
    ht->count++;
    return pos;
  }

  // Lineær probing
  do {
    ht->collisions++;
    pos = (pos + 1) % ht->size;
    if (ht->table[pos] == -1) {
      ht->table[pos] = key;
      ht->count++;
      return pos;
    }
  } while (pos != start_pos);

  return -1; // ingen ledig
}

int insert_double(HashTableV2 *ht, int key) {
  unsigned int h1 = hash1(key, ht->size);
  unsigned int pos = h1;

  if (ht->table[pos] == -1) {
    ht->table[pos] = key;
    ht->count++;
    return pos;
  }

  unsigned int h2 = hash2(key, ht->size);
  unsigned int start_pos = pos;

  if (h2 == 0) h2 = 1;

  do {
    ht->collisions++;
    pos = (pos + h2) % ht->size;

    if (ht->table[pos] == -1) {
      ht->table[pos] = key;
      ht->count++;
      return pos;
    }
  } while (pos != start_pos);

  return -1; // ingen ledig
}

void free_open_hash_table(HashTableV2 *ht) {
  free(ht->table);
  free(ht);
}

void calculate_time_difference(struct timespec start, struct timespec end, int *sec, long *nsec) {
  *sec = (int)(end.tv_sec - start.tv_sec);
  *nsec = end.tv_nsec - start.tv_nsec;
  if (*nsec < 0) {
    --(*sec);
    *nsec += 1000000000L;
  }
}

// src: https://www.geeksforgeeks.org/dsa/program-to-find-the-next-prime-number/
// Function that returns true if n
// is prime else returns false
bool isPrime(int n)
{
  // Corner cases
  if (n <= 1)  return false;
  if (n <= 3)  return true;

  // This is checked so that we can skip
  // middle five numbers in below loop
  if (n%2 == 0 || n%3 == 0) return false;

  for (int i=5; i*i<=n; i=i+6)
    if (n%i == 0 || n%(i+2) == 0)
      return false;

  return true;
}

// Function to return the smallest
// prime number greater than N
int nextPrime(int N)
{

  // Base case
  if (N <= 1)
    return 2;

  int prime = N;
  bool found = false;

  // Loop continuously until isPrime returns
  // true for a number greater than n
  while (!found) {
    prime++;

    if (isPrime(prime))
      found = true;
  }

  return prime;
}
void del2() {
  printf("\n========== DEL 2: HASHTABELLER OG YTELSE ==========\n\n");
  int N = 13000000;
  int M = nextPrime(N);

  printf("Genererer %d unike tilfeldige tall...\n", N);

  // unike tall i array
  int *numbers = malloc(N * sizeof(int));
  numbers[0] = rand() % 1000;
  for (int i = 1; i < N; i++) {
    numbers[i] = numbers[i-1] + 1 + rand() % 1000;
  }

  for (int i = 0; i < N; i++) {
    int j = rand() % N;
    int temp = numbers[i];
    numbers[i] = numbers[j];
    numbers[j] = temp;
  }

  printf("Tabellstørrelse: %d\n\n", M);

  double fill_rates[] = {0.3, 0.5, 0.8, 0.90, 0.99, 1.0};

  printf("%-15s %-20s %-15s %-15s %-15s\n", "Fyllingsgrad", "Type", "Kollisjoner", "Tid", "Kollisjoner/innsetting");
  printf("---------------------------------------------------------------------------------------------\n");

  for (int i = 0; i < 6; i++) {
    int n = (int)(N * fill_rates[i]);
    char tid_buf[64];

    // Linear probing
    HashTableV2 *ht_linear = create_open_hash_table(M);
    struct timespec start, end;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);

    for (int j = 0; j < n; j++) {
      insert_linear(ht_linear, numbers[j]);
    }

    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end);
    int sec;
    long nsec;
    calculate_time_difference(start, end, &sec, &nsec);
    format_tid(tid_buf, sizeof(tid_buf), sec, nsec);

    printf("%-15.0f%% %-20s %-15llu %-15s %-15.2f\n",
           fill_rates[i] * 100, "Lineær probing", ht_linear->collisions, tid_buf, (double) ht_linear->collisions / ht_linear->count);
    free_open_hash_table(ht_linear);

    // Double hashing
    HashTableV2 *ht_double = create_open_hash_table(M);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
    for (int j = 0; j < n; j++) {
      insert_double(ht_double, numbers[j]);
    }
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end);
    calculate_time_difference(start, end, &sec, &nsec);
    format_tid(tid_buf, sizeof(tid_buf), sec, nsec);
    printf("%-15s %-20s %-15llu %-15s %-15.2f\n",
           "", "Dobbel hashing", ht_double->collisions, tid_buf, (double) ht_double->collisions / ht_double->count);
    free_open_hash_table(ht_double);
    printf("\n");
  }

  free(numbers);
}

int main(void) {
  srand(time(NULL));
  del1();
  del2();
  return 0;
}