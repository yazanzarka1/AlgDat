#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define COMPRESS 0
#define DECOMPRESS 1
#define BLOCKSIZE (900 * 1024)  // 900 KB


typedef struct Suffix {
  int32_t index;
  int32_t rank[2];
} suffix_t;

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
  unsigned char mtf[256];
  for (int i = 0; i < 256; i++) {
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
int32_t btw_transform(const unsigned char *input_buffer,
                      const size_t input_size,
                      unsigned char *output_buffer) {
  if (input_size == 0) return -1;
  suffix_t *suffixes = malloc(input_size * sizeof(suffix_t));

  for (int32_t i = 0; i < input_size; i++) {
    suffixes[i].index = i;
    suffixes[i].rank[0] = input_buffer[i];
    suffixes[i].rank[1] = (i + 1 < input_size) ? input_buffer[i + 1] : -1;
  }

  qsort(suffixes, input_size, sizeof(suffix_t), cpm_suffix);

  for (size_t k = 4; k < 2 * input_size; k *= 2) {
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


    for (size_t i = 1; i < input_size; i++) {
      size_t next_idx = suffixes[i].index + k;
      suffixes[i].rank[0] = temp_rank[suffixes[i].index];
      suffixes[i].rank[1] = (next_idx < input_size) ? temp_rank[next_idx] : -1;
    }

    free(temp_rank);
    qsort(suffixes, input_size, sizeof(suffix_t), cpm_suffix);
  }

  int32_t original_index = -1;
  for (size_t i = 0; i < input_size; i++) {
    int32_t idx = suffixes[i].index;
    if (idx == 0) {
      original_index = i;
      output_buffer[i] = input_buffer[input_size - 1];
    } else {
      output_buffer[i] = input_buffer[idx - 1];
    }
  }
  free(suffixes);
  return original_index;
}

void reverse_btw_transform(const unsigned char *input_buffer,
                           const size_t input_size,
                           int32_t original_index,
                           unsigned char *output_buffer) {
  int32_t *T = malloc(input_size * sizeof(int32_t));
  int32_t *count = malloc(256 * sizeof(int32_t));

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


  // her mapper vi forekomster i btw-bufferen til tilsvarende forekomster i den sorterte rotasjons-arrayen
  for (size_t i = 0; i < input_size; i++) {
    T[count[input_buffer[i]]] = (int32_t) i;
    count[input_buffer[i]]++;
  }

  int32_t idx = original_index;
  for (size_t i = 0; i < input_size; i++) {
    output_buffer[i] = input_buffer[idx];
    idx = T[idx];
  }

  free(T);
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
    if (buffer[i] == buffer[i + 1]) {
      int16_t count = 1;
      while (i < buffer_size && buffer[i + count] == c && count < INT16_MAX) {
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
      while (i < buffer_size && buffer[i] != buffer[i + 1] && count > INT16_MIN) {
        count--;
        i++;
      }
      output_buffer[output_pos++] = count & 0xFF;
      output_buffer[output_pos++] = (count >> 8) & 0xFF;
      for (int32_t j = 0; j < -count; j++) {
        output_buffer[output_pos++] = buffer[i + j - (-count)];
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
  while (i < buffer_size) {
    int16_t count = (int16_t) (buffer[i] | buffer[i + 1] << 8);
    i += 2;
    if (count > 0) {
      // Repeterende sekvens
      const unsigned char c = buffer[i++];
      for (int16_t j = 0; j < count; j++) {
        output_buffer[output_pos++] = c;
      }
    } else {
      // Ikke-repeterende sekvens
      for (int16_t j = 0; j < -count; j++) {
        output_buffer[output_pos++] = buffer[i++];
      }
    }
  }
  return output_pos;
}


size_t run_compress_sequence(const unsigned char *input_buffer,
                             const size_t input_size,
                             unsigned char *output_buffer) {
  unsigned char *rle_buffer = malloc(BLOCKSIZE * 2);
  unsigned char *btw_buffer = malloc(BLOCKSIZE);
  unsigned char *mtf_buffer = malloc(BLOCKSIZE);
  size_t rle_size = run_length_encoding(input_buffer, input_size, rle_buffer);
  int32_t original_index = btw_transform(rle_buffer, rle_size, btw_buffer);
  move_to_front_encode(btw_buffer, rle_size, mtf_buffer);
  rle_size = run_length_encoding(mtf_buffer, rle_size, rle_buffer);
  // Først lagrer vi original_index som 4 bytes i output_buffer
  output_buffer[0] = (original_index >> 24) & 0xFF;
  output_buffer[1] = (original_index >> 16) & 0xFF;
  output_buffer[2] = (original_index >> 8) & 0xFF;
  output_buffer[3] = original_index & 0xFF;
  memcpy(output_buffer + 4, rle_buffer, rle_size);
  free(rle_buffer);
  free(btw_buffer);
  free(mtf_buffer);
  return rle_size + 4;
}

size_t run_decompress_sequence(const unsigned char *input_buffer,
                               const size_t input_size,
                               unsigned char *output_buffer) {
  unsigned char *rle_buffer = malloc(BLOCKSIZE * 2);
  unsigned char *btw_buffer = malloc(BLOCKSIZE);
  unsigned char *mtf_buffer = malloc(BLOCKSIZE);
  // Hent original_index fra de første 4 bytene i input_buffer
  int32_t original_index = (input_buffer[0] << 24) |
                           (input_buffer[1] << 16) |
                           (input_buffer[2] << 8) |
                           input_buffer[3];
  size_t rle_size = run_length_decoding(input_buffer + 4, input_size - 4, rle_buffer);
  move_to_front_decode(rle_buffer, rle_size, mtf_buffer);
  reverse_btw_transform(mtf_buffer, rle_size, original_index, btw_buffer);
  size_t output_size = run_length_decoding(btw_buffer, rle_size, output_buffer);
  free(rle_buffer);
  free(btw_buffer);
  free(mtf_buffer);
  return output_size;
}

size_t run_bzip2_compression(FILE *input_file,
                             FILE *output_file) {
  unsigned char *buffer = malloc(BLOCKSIZE);
  unsigned char *output_buffer = malloc(BLOCKSIZE * 2);
  size_t total_output_size = 0;
  int num_blocks = 0;
  while (!feof(input_file)) {
    ++num_blocks;
    const size_t bytes_read = fread(buffer, 1, BLOCKSIZE, input_file);
    if (bytes_read == 0) {
      break;
    }
    const size_t compressed = run_compress_sequence(buffer, bytes_read, output_buffer);
    total_output_size += compressed;
    fwrite(output_buffer, 1, compressed, output_file);
  }

  free(buffer);
  free(output_buffer);
  return total_output_size;
}

size_t run_bzip2_decompression(FILE *input_file,
                               FILE *output_file) {
  unsigned char *buffer = malloc(BLOCKSIZE);
  unsigned char *output_buffer = malloc(BLOCKSIZE * 2);
  size_t total_output_size = 0;
  int num_blocks = 0;
  while (!feof(input_file)) {
    ++num_blocks;
    const size_t bytes_read = fread(buffer, 1, BLOCKSIZE, input_file);
    if (bytes_read == 0) {
      break;
    }
    const size_t compressed = run_decompress_sequence(buffer, bytes_read, output_buffer);
    total_output_size += compressed;
    fwrite(output_buffer, 1, compressed, output_file);
  }

  free(buffer);
  free(output_buffer);
  return total_output_size;
}


int main(int argc, char *argv[]) {
  if (argc != 4) {
    print_instructions();
    return 1;
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


  if (mode == COMPRESS) {
    run_bzip2_compression(input_file, output_file);
  } else if (mode == DECOMPRESS) {
    run_bzip2_decompression(input_file, output_file);
  }


  return 0;
}
