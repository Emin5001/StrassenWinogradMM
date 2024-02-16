#include <stdlib.h> 
#include <stdio.h>

void print_matrix(double *a, int size) {
  for (int i = 0; i < size * size; i++) {
    if (i % size == 0) {
       printf("\n");
    }
    printf("%d ", (int) a[i]);
  }
}

void multiply(double* a, double* b, double* res, int size) {
  double sum = 0.0;
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      sum = 0.0;
      for (int k = 0; k < size; k++) {
        sum += a[i * size + k] * b[k * size + j];
//        printf("sum += a[%d] * b[%d]\n", i * size + k, k * size + j);
      }
      res[i * size + j] += sum;
 //     printf("c[%d] += sum\n", i * size + j);
    }
  }
}

int main() {
  double *a = (double *) malloc(9 * 9 * sizeof(double));
  double *b = (double *) malloc(9 * 9 * sizeof(double));
  for (int i = 0; i < 9 * 9; i++) {
    a[i] = i;
    b[i] = i;
  }
  // pad it to 16.
  double *padded_a = (double *) calloc (16 * 16, sizeof(double));
  double *padded_b = (double *) calloc (16 * 16, sizeof(double));
  for (int i = 0; i < 9 * 9; i++) {
    padded_a[i] = a[i];
    padded_b[i] = b[i];
  }
  printf("\npadded_aa is:\n" );
  print_matrix(padded_a, 16);
  printf("\npadded_b is:\n");
  print_matrix(padded_b, 16);
  double *res = (double *) calloc (16 * 16, sizeof(double));
  multiply(padded_a, padded_b, res, 16);
  printf("\nresulting matrix is: \n");
  print_matrix(res, 16);
  return 0;
}
