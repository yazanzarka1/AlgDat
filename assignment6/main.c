#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define COMPRESS 0
#define DECOMPRESS 1
#define BLOCKSIZE (900 * 1024)
#define ALPHABET_SIZE 256


typedef struct Suffix {
  int32_t index;
  int32_t rank[2];
} suffix_t;

typedef struct BlockHeader {
  uint32_t original_index;
  uint32_t block_size;
  uint32_t rle_size;
  uint32_t freq_table[ALPHABET_SIZE];
} block_header_t;


typedef struct PriorityQueueNode {
  int16_t byte;
  uint32_t freq;
  struct PriorityQueueNode *left, *right;
} pq_node_t;

typedef struct PriorityQueue {
  pq_node_t **data;
  size_t size;
  size_t capacity;
} pq_t;

typedef struct {
  uint64_t bits;
  uint8_t nbits;
} code_t;

typedef struct {
  code_t code[ALPHABET_SIZE];
} code_table_t;

struct BitWriterOrReader {
  uint8_t bit_buffer;
  uint8_t bit_count;
  unsigned char *buffer;
};

typedef struct BitWriterOrReader bit_writer_t;
typedef struct BitWriterOrReader bit_reader_t;

void init_bit_writer_or_reader(bit_writer_t *bw, unsigned char *outbuffer) {
  bw->bit_buffer = 0;
  bw->bit_count = 0;
  bw->buffer = outbuffer;
}

void br_put_bit(bit_writer_t *bw, uint64_t bits, uint8_t nbits) {
  for (int i = nbits - 1; i >= 0; i--) {
    bw->bit_buffer = (bw->bit_buffer << 1) | ((bits >> i) & 1);
    bw->bit_count++;
    if (bw->bit_count == 8) {
      *(bw->buffer++) = bw->bit_buffer;
      bw->bit_buffer = 0;
      bw->bit_count = 0;
    }
  }
}

void br_get_bit(bit_reader_t *br, uint8_t *bit) {
  if (br->bit_count == 0) {
    br->bit_buffer = *(br->buffer++);
    br->bit_count = 8;
  }
  *bit = (br->bit_buffer >> (br->bit_count - 1)) & 1;
  br->bit_count--;
}

void br_flush(bit_writer_t *bw) {
  if (bw->bit_count > 0) {
    bw->bit_buffer <<= (8 - bw->bit_count);
    *(bw->buffer++) = bw->bit_buffer;
    bw->bit_buffer = 0;
    bw->bit_count = 0;
  }
}


// Heap og prioritetskø implementasjon for Huffman-koding
void pq_swap_node(pq_node_t **a, pq_node_t **b) {
  pq_node_t *t = *a;
  *a = *b;
  *b = t;
}
pq_node_t *create_pq_node(int16_t byte, uint32_t freq, pq_node_t *left, pq_node_t *right) {
  pq_node_t *node = (pq_node_t *) malloc(sizeof(pq_node_t));
  node->byte = byte;
  node->freq = freq;
  node->left = left;
  node->right = right;
  return node;
}

pq_t *create_pq(size_t capacity) {
  pq_t *pq = (pq_t *) malloc(sizeof(pq_t));
  pq->data = (pq_node_t **) malloc(capacity * sizeof(pq_node_t *));
  pq->size = 0;
  pq->capacity = capacity;
  return pq;
}

void heapify(pq_t *pq, size_t index) {
  while (index > 0) {
    size_t parent = (index - 1) / 2;
    if (pq->data[parent]->freq <= pq->data[index]->freq) {
      break;
    }
    pq_swap_node(&pq->data[parent], &pq->data[index]);
    index = parent;
  }
}

void free_pq(pq_t *pq) {
  for (size_t i = 0; i < pq->size; i++) {
    free(pq->data[i]);
  }
  free(pq->data);
  free(pq);
}


void pq_push(pq_t *pq, pq_node_t *node) {
  if (pq->size >= pq->capacity) {
    pq->capacity *= 2;
    pq_node_t **temp = pq->data;
    pq->data = (pq_node_t **) realloc(pq->data, pq->capacity * sizeof(pq_node_t *));
    if (pq->data == NULL) {
      free(temp);
      fprintf(stderr, "Memory reallocation failed\n");
      exit(EXIT_FAILURE);
    }
  }
  pq->data[pq->size] = node;
  heapify(pq, pq->size++);
}

pq_node_t *pq_pop(pq_t *pq) {
  if (pq->size == 0) {
    return NULL;
  }
  pq_node_t *min_node = pq->data[0];
  pq->data[0] = pq->data[--pq->size];
  size_t index = 0;
  while (1) {
    size_t left = 2 * index + 1;
    size_t right = 2 * index + 2;
    size_t smallest = index;
    if (left < pq->size && pq->data[left]->freq < pq->data[smallest]->freq) {
      smallest = left;
    }
    if (right < pq->size && pq->data[right]->freq < pq->data[smallest]->freq) {
      smallest = right;
    }
    if (smallest == index) {
      break;
    }
    pq_swap_node(&pq->data[index], &pq->data[smallest]);
    index = smallest;
  }
  return min_node;
}


// Bygger Huffman-treet basert på frekvenstabellen

uint32_t *build_frequency_table(const unsigned char *data, size_t size) {
  uint32_t *freq_table = calloc(ALPHABET_SIZE, sizeof(uint32_t));
  for (size_t i = 0; i < size; i++) {
    freq_table[data[i]]++;
  }
  return freq_table;
}

pq_node_t *build_huffman_tree(const uint32_t *freq_table) {
  pq_t *pq = create_pq(ALPHABET_SIZE);
  for (int i = 0; i < ALPHABET_SIZE; i++) {
    if (freq_table[i] > 0) {
      pq_push(pq, create_pq_node((unsigned char) i, freq_table[i], NULL, NULL));
    }
  }

  while (pq->size > 1) {
    pq_node_t *left = pq_pop(pq);
    pq_node_t *right = pq_pop(pq);
    pq_node_t *internal = create_pq_node(-1, left->freq + right->freq, left, right);
    internal->left = left;
    internal->right = right;
    pq_push(pq, internal);
  }

  pq_node_t *root = pq_pop(pq);
  free_pq(pq);
  return root;
}


void free_tree(pq_node_t *n) {
  if (!n) return;
  free_tree(n->left);
  free_tree(n->right);
  free(n);
}

void build_code_table(const pq_node_t *node,
                      code_table_t *code_table,
                      uint64_t code,
                      uint8_t depth) {
  ;
  if (node->left == NULL && node->right == NULL) {
    code_table->code[(unsigned char) node->byte].bits = code;
    // antall biter i koden tilsvarer dybden i treet
    code_table->code[(unsigned char) node->byte].nbits = depth;
    return;
  }
  if (node->left) {
    // legge til en bit med verdi 0 for venstre barn
    build_code_table(node->left, code_table, (code << 1), depth + 1);
  }
  if (node->right) {
    // legg til en bit med verdi 1 for høyre barn
    build_code_table(node->right, code_table, (code << 1) | 1, depth + 1);
  }
}


// huffman koding og dekoding funksjoner

size_t huffman_encode(const unsigned char *input_buffer,
                      const size_t input_size,
                      const code_table_t *code_table,
                      unsigned char *output_buffer) {
  bit_writer_t bw;
  init_bit_writer_or_reader(&bw, output_buffer);
  for (size_t i = 0; i < input_size; i++) {
    unsigned char byte = input_buffer[i];
    code_t code = code_table->code[byte];
    br_put_bit(&bw, code.bits, code.nbits);
  }
  br_flush(&bw);
  return (size_t) (bw.buffer - output_buffer);
}


size_t huffman_decode(const unsigned char *input_buffer,
                      const size_t input_size,
                      const pq_node_t *huffman_tree,
                      unsigned char *output_buffer,
                      size_t output_size
) {
  bit_reader_t br;
  init_bit_writer_or_reader(&br, (unsigned char *) input_buffer);
  size_t out_pos = 0;
  const pq_node_t *current = huffman_tree;
  for (size_t i = 0; i < input_size * 8; i++) {
    uint8_t bit;
    br_get_bit(&br, &bit);
    if (bit == 0) {
      current = current->left;
    } else {
      current = current->right;
    }

    if (current->left == NULL && current->right == NULL) {
      output_buffer[out_pos++] = (unsigned char) current->byte;
      current = huffman_tree;
      if (out_pos >= output_size) break;
    }
  }
  return out_pos;
}

void print_instructions() {
  printf("Usage: main6 [command] [input_file] [output_file]\n");
  printf("Commands:\n");
  printf("  -c : Compress the input file\n");
  printf("  -d : Decompress the input file\n");
  printf("Example:\n");
  printf("  main6 -c input.txt output.bz2  # Compress input.txt to output.bz2\n");
  printf("  main6 -d output.bz2 input.txt  # Decompress output.bz2 to input.txt\n");
}


// Move-to-Front encoding
void move_to_front_encode(unsigned char *data, size_t size, unsigned char *output) {
  unsigned char mtf[ALPHABET_SIZE];
  for (int i = 0; i < ALPHABET_SIZE; i++) {
    mtf[i] = i;
  }

  for (size_t i = 0; i < size; i++) {
    unsigned char c = data[i];

    // Finn posisjon i MTF liste
    int pos = 0;
    while (mtf[pos] != c) pos++;

    output[i] = pos;

    // Flytt til front
    unsigned char temp = mtf[pos];
    memmove(&mtf[1], &mtf[0], pos);
    mtf[0] = temp;
  }
}

void move_to_front_decode(unsigned char *data, size_t size, unsigned char *output) {
  unsigned char mtf[256];

  for (int i = 0; i < 256; i++) {
    mtf[i] = i;
  }

  for (size_t i = 0; i < size; i++) {
    int pos = data[i];
    unsigned char c = mtf[pos];
    output[i] = c;

    memmove(&mtf[1], &mtf[0], pos);
    mtf[0] = c;
  }
}

int cpm_suffix(const void *a, const void *b) {
  suffix_t *suffixA = (suffix_t *) a;
  suffix_t *suffixB = (suffix_t *) b;
  if (suffixA->rank[0] == suffixB->rank[0]) {
    return suffixA->rank[1] - suffixB->rank[1];
  }
  return suffixA->rank[0] - suffixB->rank[0];
}

// ref: https://www.geeksforgeeks.org/dsa/suffix-array-set-2-a-nlognlogn-algorithm/
// Måtte referere til eksterne kilder for denne implementasjonen
// O(n log n log n)
int32_t bwt_transform(const unsigned char *input_buffer,
                      const size_t input_size,
                      unsigned char *output_buffer) {
  if (input_size == 0) return -1;
  suffix_t *suffixes = malloc(input_size * sizeof(suffix_t));

  for (int32_t i = 0; i < input_size; i++) {
    suffixes[i].index = i;
    suffixes[i].rank[0] = input_buffer[i];
    suffixes[i].rank[1] = input_buffer[(i + 1) % input_size];
  }

  qsort(suffixes, input_size, sizeof(suffix_t), cpm_suffix);

  for (size_t k = 2; k < 2 * input_size; k *= 2) {
    int32_t *temp_rank = malloc(input_size * sizeof(int32_t));
    temp_rank[suffixes[0].index] = 0;

    for (size_t i = 1; i < input_size; i++) {
      suffix_t *prev = &suffixes[i - 1];
      suffix_t *curr = &suffixes[i];
      if (curr->rank[0] == prev->rank[0] && curr->rank[1] == prev->rank[1]) {
        temp_rank[curr->index] = temp_rank[prev->index];
      } else {
        temp_rank[curr->index] = temp_rank[prev->index] + 1;
      }
    }

    for (size_t i = 0; i < input_size; i++) {
      size_t next_idx = (suffixes[i].index + k) % input_size;
      suffixes[i].rank[0] = temp_rank[suffixes[i].index];
      suffixes[i].rank[1] = temp_rank[next_idx];
    }

    free(temp_rank);
    qsort(suffixes, input_size, sizeof(suffix_t), cpm_suffix);
  }

  int32_t original_index = -1;
  for (size_t i = 0; i < input_size; i++) {
    int32_t idx = suffixes[i].index;
    if (idx == 0) {
      original_index = i;
    }
    output_buffer[i] = input_buffer[(idx + input_size - 1) % input_size];
  }
  free(suffixes);
  return original_index;
}

void reverse_bwt_transform(const unsigned char *input_buffer,
                           const size_t input_size,
                           const int32_t original_index,
                           unsigned char *output_buffer) {
  int32_t *T = malloc(input_size * sizeof(int32_t));
  unsigned char *F = malloc(input_size * sizeof(unsigned char));
  int32_t *count = malloc(256 * sizeof(int32_t));

  memset(count, 0, 256 * sizeof(int32_t));

  // her teller vi antall forekomster av en byte-verdi i input_buffer
  for (size_t i = 0; i < input_size; i++) {
    count[input_buffer[i]]++;
  }

  // her berenge vi startindeksen for hver byte-verdi i den sorterte
  // versjonen av input_buffer
  int32_t total_sum = 0;
  for (int32_t i = 0; i < 256; i++) {
    int32_t old_count = count[i];
    count[i] = total_sum;
    total_sum += old_count;
  }

  int32_t *temp_count = malloc(256 * sizeof(int32_t));
  memcpy(temp_count, count, 256 * sizeof(int32_t));

  // Bygg F (første kolonne - sortert versjon av L)
  for (size_t i = 0; i < input_size; i++) {
    int pos = temp_count[input_buffer[i]];
    F[pos] = input_buffer[i];
    T[pos] = (int32_t) i;
    temp_count[input_buffer[i]]++;
  }

  free(temp_count);

  // Rekonstruer original streng
  int32_t idx = original_index;
  for (size_t i = 0; i < input_size; i++) {
    output_buffer[i] = F[idx];
    idx = T[idx];
  }

  free(T);
  free(F);
  free(count);
}

// Run-Length Coding
size_t run_length_encoding(const unsigned char *buffer,
                           const size_t buffer_size,
                           unsigned char *output_buffer) {
  size_t output_pos = 0;
  size_t i = 0;

  while (i < buffer_size) {
    const unsigned char c = buffer[i];

    // Dersom ord[i] og ord[i+1] er like, så er det en repeterende sekvens
    // vi teller antall repeterende tegn og lagrer det i output_buffer
    if (i + 1 < buffer_size && buffer[i] == buffer[i + 1]) {
      int16_t count = 1;
      while (i + count < buffer_size && buffer[i + count] == c && count < INT16_MAX) {
        count++;
      }
      i += count;
      output_buffer[output_pos++] = count & 0xFF;
      output_buffer[output_pos++] = (count >> 8) & 0xFF;
      output_buffer[output_pos++] = c;
    } else {
      // Ikke-repeterende sekvens
      // dersom ord[i] og ord[i+1] ikke er like, teller vi antall
      // ikke-repeterende tegn og lagrer det i output_buffer med negativt fortegn
      int16_t count = 0;
      size_t start = i;
      while (i < buffer_size && (i + 1 >= buffer_size || buffer[i] != buffer[i + 1]) && count >
        INT16_MIN) {
        count--;
        i++;
      }
      output_buffer[output_pos++] = count & 0xFF;
      output_buffer[output_pos++] = (count >> 8) & 0xFF;
      for (int32_t j = 0; j < -count; j++) {
        output_buffer[output_pos++] = buffer[start + j];
      }
    }
  }

  return output_pos;
}

size_t run_length_decoding(const unsigned char *buffer,
                           const size_t buffer_size,
                           unsigned char *output_buffer) {
  size_t output_pos = 0;
  size_t i = 0;
  while (i + 1 < buffer_size) {
    int16_t count = (int16_t) (buffer[i] | buffer[i + 1] << 8);
    i += 2;
    if (count > 0) {
      // Repeterende sekvens
      if (i >= buffer_size) break;
      const unsigned char c = buffer[i++];
      for (int16_t j = 0; j < count; j++) {
        output_buffer[output_pos++] = c;
      }
    } else {
      // Ikke-repeterende sekvens
      for (int j = 0; j < -count; j++) {
        if (i >= buffer_size) break;
        output_buffer[output_pos++] = buffer[i++];
      }
    }
  }
  return output_pos;
}

void print_huffman_tree(pq_node_t *node, int depth) {
  if (node->left == NULL && node->right == NULL) {
    printf("Byte: %d, Depth: %d, Freq: %d\n", node->byte, depth, node->freq);
    return;
  }

  if (node->left) {
    print_huffman_tree(node->left, depth + 1);
  }
  if (node->right) {
    print_huffman_tree(node->right, depth + 1);
  }
}


size_t run_compress_sequence(const unsigned char *input_buffer,
                             const size_t input_size,
                             unsigned char *output_buffer) {
  unsigned char *rle_buffer = malloc(BLOCKSIZE * 2);
  unsigned char *bwt_buffer = malloc(BLOCKSIZE * 2);
  unsigned char *mtf_buffer = malloc(BLOCKSIZE * 2);
  unsigned char *rle2_buffer = malloc(BLOCKSIZE * 2);
  unsigned char *huffman_buffer = malloc(BLOCKSIZE * 2);


  size_t rle_size = run_length_encoding(input_buffer, input_size, rle_buffer);
  int32_t original_index = bwt_transform(rle_buffer, rle_size, bwt_buffer);
  move_to_front_encode(bwt_buffer, rle_size, mtf_buffer);
  size_t rle2_size = run_length_encoding(mtf_buffer, rle_size, rle2_buffer);
  // Bygg frekvenstabell for Huffman-koding
  uint32_t *freq_table = build_frequency_table(rle2_buffer, rle2_size);
  pq_node_t *huffman_tree = build_huffman_tree(freq_table);
  code_table_t *code_table = malloc(sizeof(code_table_t));
  build_code_table(huffman_tree, code_table, 0, 0);
  size_t final_size = huffman_encode(rle2_buffer,
                                     rle2_size,
                                     code_table,
                                     huffman_buffer);

  // huffman
  block_header_t header = {
    .original_index = original_index,
    .block_size = final_size,
    .rle_size = rle2_size
  };

  memcpy(header.freq_table, freq_table, ALPHABET_SIZE * sizeof(uint32_t));
  memcpy(output_buffer, &header, sizeof(block_header_t));
  memcpy(output_buffer + sizeof(block_header_t), huffman_buffer, final_size);

  free(rle_buffer);
  free(bwt_buffer);
  free(mtf_buffer);
  free(rle2_buffer);
  free(huffman_buffer);
  return final_size + sizeof(block_header_t);
}

size_t run_decompress_sequence(
  const block_header_t *header,
  const unsigned char *input_buffer,
  const size_t input_size,
  unsigned char *output_buffer) {
  pq_node_t *huffman_tree = build_huffman_tree(header->freq_table);
  unsigned char *huffman_decoded = malloc(BLOCKSIZE);
  unsigned char *rle2_decoded = malloc(BLOCKSIZE * 2);
  unsigned char *mtf_decoded = malloc(BLOCKSIZE * 2);
  unsigned char *bwt_decoded = malloc(BLOCKSIZE * 2);

  size_t huffman_size = huffman_decode(input_buffer,
                                       input_size,
                                       huffman_tree,
                                       huffman_decoded,
                                       header->rle_size);

  size_t mtf_size = run_length_decoding(huffman_decoded, huffman_size, rle2_decoded);
  move_to_front_decode(rle2_decoded, mtf_size, mtf_decoded);
  reverse_bwt_transform(mtf_decoded,
                        mtf_size,
                        (int32_t)header->original_index,
                        bwt_decoded);

  size_t output_size = run_length_decoding(bwt_decoded, mtf_size, output_buffer);


  free_tree(huffman_tree);
  free(huffman_decoded);
  free(rle2_decoded);
  free(mtf_decoded);
  free(bwt_decoded);

  return output_size;
}

size_t run_bzip2_compression(FILE *input_file,
                             FILE *output_file) {
  unsigned char *buffer = malloc(BLOCKSIZE);
  unsigned char *output_buffer = malloc(BLOCKSIZE * 2);
  size_t total_input_size = 0;
  size_t total_output_size = 0;
  int num_blocks = 0;
  fseek(input_file, 0, SEEK_END);
  printf("Estimated blocks count : %zu\n", (size_t) (ftell(input_file) / BLOCKSIZE) + 1);

  fseek(input_file, 0, SEEK_SET);
  while (!feof(input_file)) {
    const size_t bytes_read = fread(buffer, 1, BLOCKSIZE, input_file);

    if (bytes_read == 0) {
      break;
    }
    ++num_blocks;
    total_input_size += bytes_read;
    if (buffer == NULL) {
      break;
    }
    const size_t compressed = run_compress_sequence(buffer, bytes_read, output_buffer);
    total_output_size += compressed;
    fwrite(output_buffer, 1, compressed, output_file);
  }


  printf("\n=== Compression Statistics ===\n");
  printf("Blocks processed: %d\n", num_blocks);
  printf("Input size:       %zu bytes (%.2f KB)\n",
         total_input_size,
         (double) total_input_size / 1024.0);
  printf("Output size:      %zu bytes (%.2f KB)\n",
         total_output_size,
         (double) total_output_size / 1024.0);
  printf("Compression ratio: %.2f%%\n",
         (1.0 - (double) total_output_size / (double) total_input_size) * 100);

  free(buffer);
  free(output_buffer);
  return total_output_size;
}

size_t run_bzip2_decompression(FILE *input_file,
                               FILE *output_file) {
  block_header_t *block_header_buffer = malloc(sizeof(block_header_t));
  unsigned char *buffer = malloc(BLOCKSIZE * 2);
  unsigned char *output_buffer = malloc(BLOCKSIZE * 6);
  size_t total_input_size = 0;
  size_t total_output_size = 0;
  int num_blocks = 0;
  while (!feof(input_file)) {
    fread(block_header_buffer, 1, sizeof(block_header_t), input_file);
    const size_t bytes_read = fread(buffer, 1, block_header_buffer->block_size, input_file);
    if (bytes_read == 0) {
      break;
    }
    ++num_blocks;
    total_input_size += bytes_read;
    if (buffer == NULL) {
      break;
    }
    const size_t decompressed = run_decompress_sequence(block_header_buffer,
                                                        buffer,
                                                        bytes_read,
                                                        output_buffer);
    total_output_size += decompressed;
    fwrite(output_buffer, 1, decompressed, output_file);
  }

  printf("\n=== Decompression Statistics ===\n");
  printf("Blocks processed: %d\n", num_blocks);
  printf("Input size:       %zu bytes (%.2f KB)\n",
         total_input_size,
         (double) total_input_size / 1024.0);
  printf("Output size:      %zu bytes (%.2f KB)\n",
         total_output_size,
         (double) total_output_size / 1024.0);
  printf("Expansion ratio:  %.2f%%\n",
         ((double) total_output_size / (double) total_input_size - 1.0) * 100);

  free(buffer);
  free(block_header_buffer);
  free(output_buffer);
  return total_output_size;
}


int main(int argc, char *argv[]) {
  if (argc != 4) {
    print_instructions();
  }

  const char *command = argv[1];
  char *input_file_path = argv[2];
  char *output_file_path = argv[3];
  int mode;

  if (strcmp(command, "-c") == 0) {
    printf("Compressing %s to %s...\n", input_file_path, output_file_path);
    mode = COMPRESS;
  } else if (strcmp(command, "-d") == 0) {
    printf("Decompressing %s to %s...\n", input_file_path, output_file_path);
    mode = DECOMPRESS;
  } else {
    print_instructions();
    return 1;
  }

  FILE *input_file = fopen(input_file_path, "rb");
  if (input_file == NULL) {
    perror("Error opening input file");
    return 1;
  }

  FILE *output_file = fopen(output_file_path, "wb");
  if (output_file == NULL) {
    perror("Error opening output file");
    fclose(input_file);
    return 1;
  }

  switch (mode) {
    case COMPRESS:
      run_bzip2_compression(input_file, output_file);
      break;
    case DECOMPRESS:
      run_bzip2_decompression(input_file, output_file);
      break;
    default: ;
  }


  fclose(input_file);
  fclose(output_file);

  return 0;
}
