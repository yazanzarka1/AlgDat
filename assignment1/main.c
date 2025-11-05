#include <stdio.h>
#include <stdlib.h>
#include <time.h>

typedef struct {
    int profit;
    int buy_day;
    int sell_day;
} Result;

Result calculate_max_profit_from_changes(const int *changes, int size) {
    Result result = {0, -1, -1};

    int best_profit = 0;
    int current_profit = 0;
    int buy_day = 0;
    int sell_day = -1;
    int current_buy_day = 0;

    for (int i = 0; i < size; ++i) {
        if (current_profit <= 0) {
            current_buy_day = i;
            current_profit = changes[i];
        } else {
            current_profit += changes[i];
        }

        if (current_profit > best_profit) {
            best_profit = current_profit;
            buy_day = current_buy_day;
            sell_day = i + 1;
        }
    }

    result.profit = best_profit;
    result.buy_day = buy_day;
    result.sell_day = sell_day;
    return result;
}

int *create_random_array(int n, int low, int high) {
    int *arr = malloc(n * sizeof(int));
    srand((unsigned int) time(NULL));

    for (int i = 0; i < n; ++i) {
        arr[i] = rand() % (high - low + 1) + low;
    }
    return arr;
}

int main() {
    const int changes[] = {-1, 3, -9, 2, 2, -1, 2, -1, -5};
    const int changes_size = sizeof(changes) / sizeof(changes[0]);
    const Result book_sample = calculate_max_profit_from_changes(changes, changes_size);
    printf("Max profit: %d, Buy on day %d, Sell on day %d\n",
           book_sample.profit, book_sample.buy_day, book_sample.sell_day);

    const int n = 1000;
    for (int i = 1; i <= 100; i = i << 1) {
        double total_sec = 0.0;
        long total_nsec = 0;

        for (int j = 0; j < 10; ++j) {
            struct timespec start_time, end_time;
            int *ptr = create_random_array(i * n, -10, 10);
            clock_gettime(CLOCK_THREAD_CPUTIME_ID, &start_time);
            calculate_max_profit_from_changes(ptr, i * n);
            clock_gettime(CLOCK_THREAD_CPUTIME_ID, &end_time);
            free(ptr);

            /*  Regne ut tid. Vi får tidsforbruket i _sec og _nsec, fordi
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
        printf("Array size %d (%d*n) ", i * n, i);
        if (avg_sec >= 1.0) printf("%.3f sekunder", avg_sec + avg_nsec / 1000000000.0);
        else if (avg_nsec > 50000000) printf("%.1f ms", avg_nsec / 1000000.0);
        else if (avg_nsec > 50000) printf("%.1f µs", avg_nsec / 1000.0);
        else printf("%ld ns", avg_nsec);
        printf(" on average after %d iterations \n", 10);
    }

    return 0;
}
