#include <stdio.h>
#include <stdlib.h>

// Function to merge two halves
void merge(int arr[], int left, int mid, int right) {
  int n1 = mid - left + 1; // size of left half
  int n2 = right - mid; // size of right half

  // Temporary arrays
  int *L = (int *) malloc(n1 * sizeof(int));
  int *R = (int *) malloc(n2 * sizeof(int));

  // Copy data into temp arrays
  for (int i = 0; i < n1; i++)
    L[i] = arr[left + i];
  for (int j = 0; j < n2; j++)
    R[j] = arr[mid + 1 + j];

  // Merge the temp arrays back into arr[left..right]
  int i = 0, j = 0, k = left;

  while (i < n1 && j < n2) {
    if (L[i] <= R[j]) {
      arr[k] = L[i];
      printf("adding %d from left\n", L[i]);
      i++;
    } else {
      arr[k] = R[j];
      printf("adding %d from right\n", R[j]);
      j++;
    }
    k++;
  }

  // Copy remaining elements, if any
  while (i < n1) {
    arr[k] = L[i];
    i++;
    k++;
  }
  while (j < n2) {
    arr[k] = R[j];
    j++;
    k++;
  }

  // Free temporary arrays
  free(L);
  free(R);
}

// Recursive merge sort
void mergeSort(int arr[], int left, int right) {
  if (left < right) {
    int mid = left + (right - left) / 2;
    printf("sorting [%d - %d] [%d - %d]\n", left, mid, mid + 1, right);
    mergeSort(arr, left, mid); // Sort first half
    mergeSort(arr, mid + 1, right); // Sort second half

    merge(arr, left, mid, right); // Merge the two halves
  }
}

// Utility function to print an array
void printArray(int arr[], int size) {
  for (int i = 0; i < size; i++)
    printf("%d ", arr[i]);
  printf("\n");
}

// Main function for testing
int main() {
  int arr[] = {4, 2, 6, 1, 5, 8, 7, 4, 3};
  int size = sizeof(arr) / sizeof(arr[0]);

  printf("Original array:\n");
  printArray(arr, size);

  mergeSort(arr, 0, size - 1);

  printf("\nSorted array:\n");
  printArray(arr, size);

  return 0;
}
