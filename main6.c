#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define COMPRESS 0
#define DECOMPRESS 1
#define BLOCKSIZE (100 * 1024)
#define ALPHABET_SIZE 256


typedef struct Suffix {
  int32_t index;
  int32_t rank[2];
} suffix_t;

typedef struct BlockHeader {
  uint32_t original_index;
  uint32_t block_size;
  uint8_t reserved[8];
} block_header_t;



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

void reverse_btw_transform(const unsigned char *input_buffer,
                           const size_t input_size,
                           int32_t original_index,
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


size_t run_compress_sequence(const unsigned char *input_buffer,
                             const size_t input_size,
                             unsigned char *output_buffer) {
  unsigned char *rle_buffer = malloc(BLOCKSIZE * 2);
  unsigned char *btw_buffer = malloc(BLOCKSIZE * 2);
  unsigned char *mtf_buffer = malloc(BLOCKSIZE * 2);
  unsigned char *rle2_buffer = malloc(BLOCKSIZE * 2);

  size_t rle_size = run_length_encoding(input_buffer, input_size, rle_buffer);
  int32_t original_index = btw_transform(rle_buffer, rle_size, btw_buffer);
  move_to_front_encode(btw_buffer, rle_size, mtf_buffer);
  size_t final_size = run_length_encoding(mtf_buffer, rle_size, rle2_buffer);

  block_header_t header = {
    .original_index = original_index,
    .block_size = rle_size,
    .reserved = {0}
  };
  memcpy(output_buffer, &header, sizeof(block_header_t));
  memcpy(output_buffer + sizeof(block_header_t), rle2_buffer, final_size);

  free(rle_buffer);
  free(btw_buffer);
  free(mtf_buffer);
  free(rle2_buffer);
  return final_size + sizeof(block_header_t);
}

size_t run_decompress_sequence(const unsigned char *input_buffer,
                               const size_t input_size,
                               unsigned char *output_buffer) {
  unsigned char *rle_buffer = malloc(BLOCKSIZE * 2);
  unsigned char *btw_buffer = malloc(BLOCKSIZE * 2);
  unsigned char *mtf_buffer = malloc(BLOCKSIZE * 2);

  block_header_t header;
  memcpy(&header, input_buffer, sizeof(block_header_t));
  int32_t original_index = header.original_index;
  uint32_t bwt_size = header.block_size;

  size_t mtf_size = run_length_decoding(input_buffer + sizeof(block_header_t),
                                         input_size - sizeof(block_header_t),
                                         mtf_buffer);
  move_to_front_decode(mtf_buffer, bwt_size, btw_buffer);
  reverse_btw_transform(btw_buffer,
                        bwt_size,
                        original_index,
                        rle_buffer);
  size_t output_size = run_length_decoding(rle_buffer,
                                           bwt_size,
                                           output_buffer);

  free(rle_buffer);
  free(btw_buffer);
  free(mtf_buffer);
  return output_size;
}

size_t run_bzip2_compression(FILE *input_file,
                             FILE *output_file) {
  unsigned char *buffer = malloc(BLOCKSIZE);
  unsigned char *output_buffer = malloc(BLOCKSIZE * 2);
  size_t total_input_size = 0;
  size_t total_output_size = 0;
  int num_blocks = 0;

  clock_t start = clock();

  while (!feof(input_file)) {
    const size_t bytes_read = fread(buffer, 1, BLOCKSIZE, input_file);
    if (bytes_read == 0) {
      break;
    }
    ++num_blocks;
    total_input_size += bytes_read;

    const size_t compressed = run_compress_sequence(buffer, bytes_read, output_buffer);
    total_output_size += compressed;
    fwrite(output_buffer, 1, compressed, output_file);
  }

  clock_t end = clock();
  double elapsed = (double)(end - start) / CLOCKS_PER_SEC;

  printf("\n=== Compression Statistics ===\n");
  printf("Blocks processed: %d\n", num_blocks);
  printf("Input size:       %zu bytes (%.2f KB)\n", total_input_size, total_input_size / 1024.0);
  printf("Output size:      %zu bytes (%.2f KB)\n", total_output_size, total_output_size / 1024.0);
  printf("Compression ratio: %.2f%%\n", (1.0 - (double)total_output_size / total_input_size) * 100);
  printf("Time elapsed:     %.3f seconds\n", elapsed);
  printf("Throughput:       %.2f KB/s\n", (total_input_size / 1024.0) / elapsed);

  free(buffer);
  free(output_buffer);
  return total_output_size;
}

size_t run_bzip2_decompression(FILE *input_file,
                               FILE *output_file) {
  unsigned char *buffer = malloc(BLOCKSIZE * 2);
  unsigned char *output_buffer = malloc(BLOCKSIZE * 2);
  size_t total_input_size = 0;
  size_t total_output_size = 0;
  int num_blocks = 0;

  clock_t start = clock();

  while (!feof(input_file)) {
    const size_t bytes_read = fread(buffer, 1, BLOCKSIZE * 2, input_file);
    if (bytes_read == 0) {
      break;
    }
    ++num_blocks;
    total_input_size += bytes_read;

    const size_t decompressed = run_decompress_sequence(buffer, bytes_read, output_buffer);
    total_output_size += decompressed;
    fwrite(output_buffer, 1, decompressed, output_file);
  }

  clock_t end = clock();
  double elapsed = (double)(end - start) / CLOCKS_PER_SEC;

  printf("\n=== Decompression Statistics ===\n");
  printf("Blocks processed: %d\n", num_blocks);
  printf("Input size:       %zu bytes (%.2f KB)\n", total_input_size, total_input_size / 1024.0);
  printf("Output size:      %zu bytes (%.2f KB)\n", total_output_size, total_output_size / 1024.0);
  printf("Expansion ratio:  %.2f%%\n", ((double)total_output_size / total_input_size - 1.0) * 100);
  printf("Time elapsed:     %.3f seconds\n", elapsed);
  printf("Throughput:       %.2f KB/s\n", (total_output_size / 1024.0) / elapsed);

  free(buffer);
  free(output_buffer);
  return total_output_size;
}


int main(int argc, char *argv[]) {
  char *test_buffer =
      "Thiiiiis iiiiiiis a prepreatttijninindsf text file wiiith some errrorss.\nIt iiis meaaant to teest the codding abillityy of an AI moddel.\nPleaaase corrrect the errrorss and imprrove the quallity.";
  char *encoded_buffer = malloc(strlen(test_buffer) * 4);
  char *decoded_buffer = malloc(strlen(test_buffer) + 1);
  memset(decoded_buffer, 0, strlen(test_buffer) + 1);

  printf("=== BWT Roundtrip Test ===\n");
  size_t compressed_size = run_compress_sequence((unsigned char *)test_buffer,
                                                  strlen(test_buffer),
                                                  (unsigned char *)encoded_buffer);
  size_t decompressed_size = run_decompress_sequence((unsigned char *)encoded_buffer,
                                                      compressed_size,
                                                      (unsigned char *)decoded_buffer);

  printf("Original size:     %zu bytes\n", strlen(test_buffer));
  printf("Compressed size:   %zu bytes\n", compressed_size);
  printf("Decompressed size: %zu bytes\n", decompressed_size);
  printf("Compression ratio: %.2f%%\n\n", (1.0 - (double)compressed_size / strlen(test_buffer)) * 100);

  int match = memcmp(test_buffer, decoded_buffer, strlen(test_buffer));
  if (match == 0) {
    printf("✓ BWT roundtrip successful!\n\n");
  } else {
    printf("✗ BWT roundtrip failed!\n");
    printf("Original: %s\n", test_buffer);
    printf("Decoded:  %s\n\n", decoded_buffer);
  }

  free(encoded_buffer);
  free(decoded_buffer);

  if (argc != 4) {
    print_instructions();
    return match == 0 ? 0 : 1;
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