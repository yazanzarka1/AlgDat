#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>


long long checksum(const int *t, int n) {
  long long s = 0;
  for (int i = 0; i < n; ++i) s += (long long) t[i];
  return s;
}

void print_liste(const int *t, int n) {
  for (int i = 0; i < n; ++i) printf("%d ", t[i]);
  printf("\n");
}

int er_sortert(const int *t, int n) {
  for (int i = 1; i < n; ++i) if (t[i] < t[i - 1]) return 0;
  return 1;
}


void bytt(int *a, int *b) {
  int tmp = *a;
  *a = *b;
  *b = tmp;
}

int median3sort(int *t, int v, int h) {
  int m = (v + h) / 2;
  if (t[v] > t[m]) bytt(&t[v], &t[m]);
  if (t[m] > t[h]) {
    bytt(&t[m], &t[h]);
    if (t[v] > t[m]) bytt(&t[v], &t[m]);
  }
  return m;
}

int splitt(int *t, int v, int h) {
  int m = median3sort(t, v, h);
  int dv = t[m];
  bytt(&t[m], &t[h - 1]);
  int iv = v, ih = h - 1;
  for (;;) {
    while (t[++iv] < dv);
    while (t[--ih] > dv);
    if (iv >= ih) break;
    bytt(&t[iv], &t[ih]);
  }
  bytt(&t[iv], &t[h - 1]);
  return iv;
}


void quicksort(int *t, int v, int h) {
  if (h - v > 2) {
    int delepos = splitt(t, v, h);
    quicksort(t, v, delepos - 1);
    quicksort(t, delepos + 1, h);
  } else {
    median3sort(t, v, h);
  }
}


int plasser_min_maks(int *t, int n) {
  if (n <= 1) return 1;
  int imin = 0, imax = 0;
  int vmin = t[0], vmax = t[0];
  for (int i = 1; i < n; ++i) {
    if (t[i] < vmin) {
      vmin = t[i];
      imin = i;
    }
    if (t[i] > vmax) {
      vmax = t[i];
      imax = i;
    }
  }
  if (vmin == vmax) return 1;
  if (imin != 0) {
    bytt(&t[0], &t[imin]);
    if (imax == 0) imax = imin;
  }
  if (imax != n - 1) bytt(&t[n - 1], &t[imax]);
  return 0;
}

void forbedret_quicksort_rekursiv(int *t, int v, int h) {
  if (h - v > 2) {
    if (t[v - 1] == t[h + 1]) {
      return;
    };
    int delepos = splitt(t, v, h);
    forbedret_quicksort_rekursiv(t, v, delepos - 1);
    forbedret_quicksort_rekursiv(t, delepos + 1, h);
  } else {
    median3sort(t, v, h);
  }
}

void forbedret_quicksort(int *t, int n) {
  if (n <= 1) return;
  if (plasser_min_maks(t, n)) return;
  if (n > 2) forbedret_quicksort_rekursiv(t, 1, n - 2);
}


void fyll_mange_duplikater(int *t, int n) {
  for (int i = 0; i < n; ++i) {
    t[i] = i % 2 == 0 ? 42 : rand();
  }
}

void fyll_unike_tilfeldig(int *t, int n) {
  for (int i = 0; i < n; ++i) { t[i] = rand(); }
}

void beregn_tid_forskjell(struct timespec start, struct timespec end, int *sec, long *nsec) {
  *sec = (int) (end.tv_sec - start.tv_sec);
  *nsec = end.tv_nsec - start.tv_nsec;
  if (*nsec < 0) {
    --(*sec);
    *nsec += 1000000000L;
  }
}

void print_tid(int n, const char *suffix, int sec, long nsec) {
  if (sec >= 1.0) {
    printf("%.3f sekunder", sec + nsec / 1000000000.0);
  } else if (nsec > 50000000) {
    printf("%.1f ms", nsec / 1000000.0);
  } else if (nsec > 50000) {
    printf("%.1f µs", nsec / 1000.0);
  } else {
    printf("%ld ns", nsec);
  }
  printf("%s", suffix);
}

void valider_resultater(const int *base, int n, long long expected_checksum, const char *context) {
  if (expected_checksum != checksum(base, n)) {
    fprintf(stderr, "[%s] FEIL: sjekksum endret!\n", context);
  } else {
    printf("Checksum OK\n");
  }

  if (!er_sortert(base, n)) {
    fprintf(stderr, "[%s] FEIL: ikke sortert!\n", context);
  } else {
    printf("Sortert OK\n");
  }
}

void kjor_hash(int *base, int n, int valg, int *sec, long *nsec) {
  struct timespec st, en;

  if (valg == 1) {
    clock_gettime(CLOCK_THREAD_CPUTIME_ID, &st);
    quicksort(base, 0, n - 1);
    clock_gettime(CLOCK_THREAD_CPUTIME_ID, &en);
  }
  if (valg == 2) {
    clock_gettime(CLOCK_THREAD_CPUTIME_ID, &st);
    forbedret_quicksort(base, n);
    clock_gettime(CLOCK_THREAD_CPUTIME_ID, &en);
  }

  if (valg == 3) {
    clock_gettime(CLOCK_THREAD_CPUTIME_ID, &st);
    forbedret_quicksort(base, n);
    clock_gettime(CLOCK_THREAD_CPUTIME_ID, &en);
  }

  beregn_tid_forskjell(st, en, sec, nsec);
}

void test_metode(int valg) {
  if (valg == 1) printf("==== quicksort  ====\n");
  else if (valg == 2) printf("==== forbedret quicksort ====\n");
  else if (valg == 3) printf("==== forbedret quicksort uten duplikater ====\n");

  int sec;
  long nsec;
  // int n = 10;
  int n = 1000000;
  printf("Kjører med %d elementer \n", n);
  if (valg == 1) {
    int *base = malloc((size_t) n * sizeof(int));
    fyll_mange_duplikater(base, n);
    long long sjekk_f = checksum(base, n);
    kjor_hash(base, n, valg, &sec, &nsec);
    print_tid(n, " \n", sec, nsec);
    valider_resultater(base, n, sjekk_f, "duplikater uten forbedring");
    kjor_hash(base, n, valg, &sec, &nsec);
    print_tid(n, " sortert på nytt \n", sec, nsec);
    valider_resultater(base, n, sjekk_f, "duplikater uten forbedring igjen");
    free(base);
  }
  if (valg == 2) {
    int *base = malloc((size_t) n * sizeof(int));
    fyll_mange_duplikater(base, n);
    long long sjekk_f = checksum(base, n);
    kjor_hash(base, n, valg, &sec, &nsec);
    print_tid(n, " \n", sec, nsec);
    valider_resultater(base, n, sjekk_f, "duplikater");
    kjor_hash(base, n, valg, &sec, &nsec);
    print_tid(n, " sortert på nytt \n", sec, nsec);
    valider_resultater(base, n, sjekk_f, "duplikater igjen");
    free(base);
  }
  if (valg == 3) {
    int *base = malloc((size_t) n * sizeof(int));
    fyll_unike_tilfeldig(base, n);
    long long sjekk_f = checksum(base, n);
    kjor_hash(base, n, valg, &sec, &nsec);
    print_tid(n, " \n", sec, nsec);
    valider_resultater(base, n, sjekk_f, "unik forbedret");
    kjor_hash(base, n, valg, &sec, &nsec);
    print_tid(n, " sortert på nytt \n", sec, nsec);
    valider_resultater(base, n, sjekk_f, "unik forbedret igjen");
    free(base);
  }
}

int main(void) {
  srand((unsigned int) time(NULL));
  test_metode(1);
  test_metode(2);
  test_metode(3);
  return 0;
}
