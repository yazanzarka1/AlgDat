#include <math.h>
#include <stdio.h>
#include <time.h>
double method_1(double n, double x) {
  if (n == 1) {
    return x;
  }
  return x + method_1(n - 1, x);
}
double method_2(double n, double x) {
  if (n == 1) {
    return x;
  }
  if ((int)n % 2 == 0.0) {
    return method_2(n / 2, x + x);
  } else {
    return x + method_2((n - 1) / 2, x + x);
  }
}
void test_method(int method_choice) {

  printf("Method %d with double multiplication \n", method_choice);
  double result;
  if (method_choice == 1) {
    result = method_1(13, 2.5);
    printf("13x2.5 = %f \n", result);
    result = method_1(14, 10.1);
    printf("14x10.1 = %f\n", result);
  } else {
    result = method_2(13, 2.5);
    printf("13x2.5 = %f\n", result);
    result = method_2(14, 10.1);
    printf("14x10.1 = %f\n", result);
  }

  for (int i = 1; i < pow(2, 18); i = i << 1) {
    double total_sec = 0.0;
    long total_nsec = 0;
    for (int j = 0; j < 10; ++j) {
      struct timespec start_time, end_time;
      if (method_choice == 1) {
        clock_gettime(CLOCK_THREAD_CPUTIME_ID, &start_time);
        method_1(i, 12);
        clock_gettime(CLOCK_THREAD_CPUTIME_ID, &end_time);
      } else {
        clock_gettime(CLOCK_THREAD_CPUTIME_ID, &start_time);
        method_2(i, 12);
        clock_gettime(CLOCK_THREAD_CPUTIME_ID, &end_time);
      }
      /*  Regne ut tid. Vi får tidsforbruket i *sec og *nsec, fordi
          tid i nanosekunder blir for stort for 64-bits int.
          ref: blackboard ctiming.c
      */
      int sec = end_time.tv_sec - start_time.tv_sec;
      long nsec = end_time.tv_nsec - start_time.tv_nsec;
      if (nsec < 0) {
        sec -= 1;
        nsec += 1000000000;
      }
      total_sec += sec;
      total_nsec += nsec;
    }
    double avg_sec = total_sec / 10.0;
    long avg_nsec = total_nsec / 10;
    // ref: blackboard ctiming.c
    if (avg_nsec >= 1000000000) {
      avg_sec += avg_nsec / 1000000000;
      avg_nsec = avg_nsec % 1000000000;
    }
    /* Skriv gjennomsnittlig tidsmåling pent, uten plagsomt mange desimaler: */
    printf("%dx%d ", i, 12);
    if (avg_sec >= 1.0)
      printf("%.3f sekunder", avg_sec + avg_nsec / 1000000000.0);
    else if (avg_nsec > 50000000)
      printf("%.1f ms", avg_nsec / 1000000.0);
    else if (avg_nsec > 50000)
      printf("%.1f µs", avg_nsec / 1000.0);
    else
      printf("%ld ns", avg_nsec);
    printf(" on average after %d iterations \n", 10);
  }
}
int main() {
  printf("==== method_1 ====\n");
  test_method(1);
  printf("==== method_2 ====\n");
  test_method(2);
  return 0;
}
