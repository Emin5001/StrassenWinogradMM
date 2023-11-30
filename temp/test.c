#include <stdlib.h> 
#include <stdio.h>

void func(double* a, double* b, double* res, int size) {
  double sum = 0.0;
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      sum = 0.0;
      for (int k = 0; k < size; k++) {
        sum += a[i * size + k] * b[k * size + j];
        printf("sum += a[%d] * b[%d]\n", i * size + k, k * size + j);
      }
      res[i * size + j] += sum;
      printf("c[%d] += sum\n", i * size + j);
    }
  }
}

int main() {
  double *a = calloc(2 * 2, sizeof (double));
  double *b =  calloc(2 * 2, sizeof (double));
  double *c =  calloc(2 * 2, sizeof (double));
  func(a, b, c, 2);

  return 0;
}
