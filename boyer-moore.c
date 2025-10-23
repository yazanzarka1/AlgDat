#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define ALPHABET_SIZE 256
#define MAX(a, b) ((a) > (b) ? (a) : (b))

// ============================================================================
// BAD CHARACTER RULE
// ============================================================================
// Builds a table that shows the rightmost position (from the left)
// where each character occurs in the pattern.
// If the character is not found, the value is set to -1.
void build_bad_char_table(char *pattern, int m, int bc_table[ALPHABET_SIZE]) {
    // Initialize all to -1 (character not found)
    for (int i = 0; i < ALPHABET_SIZE; i++) {
        bc_table[i] = -1;
    }

    // Store the rightmost position for each character in the pattern
    for (int i = 0; i < m; i++) {
        bc_table[(unsigned char)pattern[i]] = i;
    }
}

// ============================================================================
// HELPER FUNCTIONS FOR GOOD SUFFIX RULE
// ============================================================================

// Checks if suffix[pos..m-1] is also a prefix of the pattern
int is_prefix(char *pattern, int m, int pos) {
    int suffix_len = m - pos;
    for (int i = 0; i < suffix_len; i++) {
        if (pattern[i] != pattern[pos + i]) {
            return 0;
        }
    }
    return 1;
}

// Finds the length of the longest suffix from position 'pos'
// that matches a suffix of the pattern
int suffix_length(char *pattern, int m, int pos) {
    int len = 0;
    int i = pos;
    int j = m - 1;

    while (i >= 0 && pattern[i] == pattern[j]) {
        len++;
        i--;
        j--;
    }
    return len;
}

// ============================================================================
// GOOD SUFFIX RULE
// ============================================================================
// Builds a table showing how far to shift the pattern
// when a mismatch occurs at position j (after matching suffix j+1..m-1)
void build_good_suffix_table(char *pattern, int m, int gs_table[]) {
    int last_prefix_pos = m;

    // Case 2: Suffix that is also a prefix
    for (int i = m - 1; i >= 0; i--) {
        if (is_prefix(pattern, m, i + 1)) {
            last_prefix_pos = i + 1;
        }
        // Distance to this prefix + remaining length
        gs_table[i] = last_prefix_pos + (m - 1 - i);
    }

    // Case 1: Suffix that occurs elsewhere in the pattern
    for (int i = 0; i < m - 1; i++) {
        int len = suffix_length(pattern, m, i);
        if (len > 0 && pattern[i - len + 1] != pattern[m - len]) {
            gs_table[m - len] = m - 1 - i + len;
        }
    }
}

// ============================================================================
// BOYER-MOORE SEARCH ALGORITHM WITH GALIL RULE
// ============================================================================
void boyer_moore_search(char *text, char *pattern) {
    int n = strlen(text);
    int m = strlen(pattern);

    if (m > n || m == 0) {
        printf("Pattern is too long or empty\n");
        return;
    }

    // Preprocess pattern
    int bc_table[ALPHABET_SIZE];
    int *gs_table = (int *)malloc(m * sizeof(int));

    build_bad_char_table(pattern, m, bc_table);
    build_good_suffix_table(pattern, m, gs_table);

    printf("BOYER-MOORE SEARCH ALGORITHM\n");

    printf("Text:     %s\n", text);
    printf("Pattern:  %s\n", pattern);

    // Show preprocessed tables
    printf("BAD CHARACTER TABLE:\n");
    for (int i = 0; i < m; i++) {
        printf("  '%c' to position %d\n", pattern[i], bc_table[(unsigned char)pattern[i]]);
    }

    printf("\nGOOD SUFFIX TABLE:\n");
    for (int i = 0; i < m; i++) {
        printf("  position %d to shift %d\n", i, gs_table[i]);
    }
    printf("\n");

    int shift = 0;
    int matches = 0;
    int galil_skip = 0;
    int iteration = 0;

    while (shift <= n - m) {
        iteration++;
        printf("ITERATION %d (shift = %d)\n", iteration, shift);
        

        // Show alignment
        printf("Text:    %s\n", text);
        printf("         ");
        for (int i = 0; i < shift; i++) printf(" ");
        printf("%s\n", pattern);

        int j = m - 1;  // Start from right end

        // Galil rule: skip known matching part
        if (galil_skip > 0 && j >= galil_skip) {
            printf("Galil: Skipping %d characters (known match)\n", galil_skip);
            j = galil_skip - 1;
        }

        printf("Comparing from right to left:\n");

        // Compare right to left
        while (j >= 0 && pattern[j] == text[shift + j]) {
            printf("  pos %d: '%c' = '%c' \n", j, pattern[j], text[shift + j]);
            j--;
        }

        if (j < 0) {
            // MATCH FOUND
            matches++;
            printf("\nMATCH #%d FOUND AT POSITION %d!\n", matches, shift);

            // Compute next shift and Galil skip
            int next_shift = gs_table[0];
            galil_skip = m - next_shift;

            printf("Using good suffix rule: shift = %d\n", next_shift);
            printf("Galil next: skip %d characters\n", galil_skip);

            shift += next_shift;
        } else {
            // MISMATCH
            printf("  pos %d: '%c' not '%c' NO MISMATCH\n",
                   j, pattern[j], text[shift + j]);

            // Compute shift based on both rules
            char mismatch_char = text[shift + j];
            int bc_shift = j - bc_table[(unsigned char)mismatch_char];
            int gs_shift = gs_table[j];

            printf("\nComputing shift:\n");
            printf("  Bad character rule:\n");
            printf("    Character '%c' in text at pattern-pos %d\n", mismatch_char, j);
            printf("    Rightmost occurrence of '%c' in pattern: %d\n",
                   mismatch_char, bc_table[(unsigned char)mismatch_char]);
            printf("    Shift = %d - %d = %d\n", j, bc_table[(unsigned char)mismatch_char], bc_shift);

            printf("  Good suffix rule:\n");
            printf("    Mismatch at pos %d, shift = %d\n", j, gs_shift);

            int max_shift = MAX(bc_shift, gs_shift);
            printf("  Choosing maximum: %d\n", max_shift);

            shift += max_shift;
            galil_skip = 0;  // Reset Galil on mismatch
        }
        printf("\n");
    }

    printf("SEARCH COMPLETE\n");
    printf("Total matches found: %d\n", matches);
    printf("Total iterations: %d\n", iteration);

    free(gs_table);
}

// ============================================================================
// EXAMPLES WITH DETAILED ANALYSIS
// ============================================================================

void example_1() {
    printf("\n\n");

    char text[] = "GCTTCTGCTACCTTTTGCGCGCGCGCGGAA";
    char pattern[] = "GCGCGC";

    printf("This example shows:\n");
    printf("- How the bad character rule works\n");
    printf("- How the pattern is shifted efficiently\n");
    printf("- Comparison from right to left\n\n");

    boyer_moore_search(text, pattern);
}

void example_2() {
    printf("\n\n");

    char text[] = "ABABABABCABABABABC";
    char pattern[] = "ABABC";

    printf("This example shows:\n");
    printf("- The Galil rule in action\n");
    printf("- Good suffix rule with overlapping patterns\n");
    printf("- Multiple matches in the same text\n\n");

    boyer_moore_search(text, pattern);
}

void example_3() {
    printf("\n\n");

    char text[] = "THIS IS A TEST TEXT FOR TESTING";
    char pattern[] = "TEST";

    printf("This example shows:\n");
    printf("- Search in natural text\n");
    printf("- Efficient skipping of irrelevant parts\n\n");

    boyer_moore_search(text, pattern);
}

// ============================================================================
// MAIN PROGRAM
// ============================================================================
int main() {
    example_1();
    example_2();
    example_3();

    return 0;
}
