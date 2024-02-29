#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>
#include <math.h>
// double** strassen (int, double **, double **, double **, double **, double **, double **, double **, double **, double **, double **, double **, double **, double **, double **, double **);

struct pad_struct {
    double **padded_a;
    double **padded_b;
    int padded_size;
    int original_size;
};

double** strassen (int, double **, double **, double **, double **, double **, 
                    double **, double **);
                    
void naive (double **, double **, double **, int);
double** matrix_add (double **, double **, double **, int);
double** matrix_subtr (double **, double **, double **, int);
void print_matrix (double **, int);
struct pad_struct pad_matrices(double **, double **, int, int);
double** unpad_matrix(double **, int, int);
bool power_of_two(int);
void jki                        (double **, double**, double**, int);
void kji                        (double **, double**, double**, int);
void ijk                        (double **, double**, double**, int);
void jik                        (double **, double**, double**, int);
void kij                        (double **, double**, double**, int);
void ikj                        (double **, double**, double**, int);

int sizes[16] = {25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 375, 400};
unsigned int TR = 2, TC = 2;

int main() {
  bool padded;
  struct pad_struct padded_matrices;
  struct timespec tick, tock;
  long long elapsed_ns;
  for (int k = 0; k < 16; k++) {
    int n = sizes[k];
    int newSize = n / 2;
    padded = false;
    double **A = (double **) malloc (n * sizeof (double *));
    double **B = (double **) malloc (n * sizeof (double *));
    double **D = (double **) malloc (n * sizeof (double *));
    for (int i = 0; i < n; i++) {
      A[i] = (double*) malloc (n * sizeof(double));
      B[i] = (double*) malloc (n * sizeof(double));
      D[i] = (double*) malloc (n * sizeof(double));
    }

    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        A[i][j] = i;
        B[i][j] = j;
      }
    }


    clock_gettime(CLOCK_REALTIME, &tick);
    if (!power_of_two(n)) { 
            double temp = (double) n / 16;
            if (temp != floor(temp)) { // if temp is NOT an integer
                TR = ceil(temp);
                TC = ceil(temp);
                padded_matrices = pad_matrices(A, B, TR, n);             
                padded = true; 
                newSize =  padded_matrices.padded_size / 2;
                n = padded_matrices.padded_size;
            } else {
                TR = temp;
                TC = temp;
            }
    } 
    newSize = n / 2;
    double **C = (double **) malloc (n * sizeof (double *));
    for (int i = 0; i < n; i++) {
        C[i] = (double *) malloc (n * sizeof(double));
    }
    double **c11 = (double **) malloc (newSize * sizeof (double *));
    double **c12 = (double **) malloc (newSize * sizeof (double *)); 
    double **c21 = (double **) malloc (newSize * sizeof (double *));
    double **c22 = (double **) malloc (newSize * sizeof (double *));

    for (int i = 0; i < n / 2; i++){
      c11[i] = (double *) malloc (newSize * sizeof (double));
      c12[i] = (double *) malloc (newSize * sizeof (double));
      c21[i] = (double *) malloc (newSize * sizeof (double));
      c22[i] = (double *) malloc (newSize * sizeof (double));
    }
    if (padded) {
        strassen (n, padded_matrices.padded_a, padded_matrices.padded_b, C, c11, c12, c21, c22);
        C = unpad_matrix(C, padded_matrices.padded_size, padded_matrices.original_size);
    } else {
        strassen(n, A, B, C, c11, c12, c21, c22);
    }
      

      free(c11);
      free(c12);
      free(c21);
      free(c22);

      clock_gettime(CLOCK_REALTIME, &tock);
      elapsed_ns = (tock.tv_sec - tick.tv_sec) * 1000000000LL + (tock.tv_nsec - tick.tv_nsec);
      printf("For size %d, Strassen elapsed_ns is %lld\n", sizes[k], elapsed_ns);
      free(A);
      free(B);    
      free(C); 
    }
    printf("Now naive timings:\n");

    for (int k = 0; k < 16; k++) {
        int n = sizes[k];
        double** A = (double **) malloc(n * sizeof(double *));
        double** B = (double **) malloc(n * sizeof(double *));
        double** C = (double **) malloc(n * sizeof(double *));

        for (int i = 0; i < n; i++) {
            A[i] = (double *) calloc (n, sizeof(double));
            B[i] = (double *) calloc (n, sizeof(double));
            C[i] = (double *) calloc (n, sizeof(double));
        }

        for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++){
                A[i][j] = i;
                B[i][j] = j;
                C[i][j] = 0;
            }
        }

        clock_gettime(CLOCK_REALTIME, &tick);
        jki(A, B, C, n);
        clock_gettime(CLOCK_REALTIME, &tock);
        elapsed_ns = (tock.tv_sec - tick.tv_sec) * 1000000000LL + (tock.tv_nsec - tick.tv_nsec);
        printf("For size %d, jki elapsed_ns is %lld\n", n, elapsed_ns);
        
        printf("***********************\n");

        clock_gettime(CLOCK_REALTIME, &tick);
        kji(A, B, C, n);
        clock_gettime(CLOCK_REALTIME, &tock);
        elapsed_ns = (tock.tv_sec - tick.tv_sec) * 1000000000LL + (tock.tv_nsec - tick.tv_nsec);
        printf("For size %d, kji elapsed_ns is %lld\n", n, elapsed_ns);

        printf("***********************\n");

        clock_gettime(CLOCK_REALTIME, &tick);
        ijk(A, B, C, n);
        clock_gettime(CLOCK_REALTIME, &tock);
        elapsed_ns = (tock.tv_sec - tick.tv_sec) * 1000000000LL + (tock.tv_nsec - tick.tv_nsec);
        printf("For size %d, ijk elapsed_ns is %lld\n", n, elapsed_ns);

        printf("***********************\n");

        clock_gettime(CLOCK_REALTIME, &tick);
        jik(A, B, C, n);
        clock_gettime(CLOCK_REALTIME, &tock);
        elapsed_ns = (tock.tv_sec - tick.tv_sec) * 1000000000LL + (tock.tv_nsec - tick.tv_nsec);
        printf("For size %d, jik elapsed_ns is %lld\n", n, elapsed_ns);

        printf("***********************\n");

        clock_gettime(CLOCK_REALTIME, &tick);
        kij(A, B, C, n);
        clock_gettime(CLOCK_REALTIME, &tock);
        elapsed_ns = (tock.tv_sec - tick.tv_sec) * 1000000000LL + (tock.tv_nsec - tick.tv_nsec);
        printf("For size %d, kij elapsed_ns is %lld\n", n, elapsed_ns);

        printf("***********************\n");

        clock_gettime(CLOCK_REALTIME, &tick);
        ikj(A, B, C, n);
        clock_gettime(CLOCK_REALTIME, &tock);
        elapsed_ns = (tock.tv_sec - tick.tv_sec) * 1000000000LL + (tock.tv_nsec - tick.tv_nsec);
        printf("For size %d, ikj elapsed_ns is %lld\n", n, elapsed_ns);
    }
    return 0;
  }

double** strassen(int n, double **A, double **B, double **C, double **c11, 
                  double **c12, double **c21, double **c22) {
    if (n == TR) {
        C[0][0] = (A[0][0] * B[0][0]) + (A[0][1] * B[1][0]);
        C[0][1] = (A[0][0] * B[0][1]) + (A[0][1] * B[1][1]);
        C[1][0] = (A[1][0] * B[0][0]) + (A[1][1] * B[1][0]);
        C[1][1] = (A[1][0] * B[0][1]) + (A[1][1] * B[1][1]);
        return C;

    }
    int newSize = n / 2;
    double **a11 = (double **) malloc (newSize * sizeof (double *));
    double **a12 = (double **) malloc (newSize * sizeof (double *));
    double **a21 = (double **) malloc (newSize * sizeof (double *));
    double **a22 = (double **) malloc (newSize * sizeof (double *));

    double **b11 = (double **) malloc (newSize * sizeof (double *));
    double **b12 = (double **) malloc (newSize * sizeof (double *));
    double **b21 = (double **) malloc (newSize * sizeof (double *));
    double **b22 = (double **) malloc (newSize * sizeof (double *));

    double **s1 = (double **) malloc (newSize * sizeof (double *));
    double **s2 = (double **) malloc (newSize * sizeof (double *));
    double **s3 = (double **) malloc (newSize * sizeof (double *));
    double **s4 = (double **) malloc (newSize * sizeof (double *));

    double **t1 = (double **) malloc (newSize * sizeof (double *));
    double **t2 = (double **) malloc (newSize * sizeof (double *));
    double **t3 = (double **) malloc (newSize * sizeof (double *));
    double **t4 = (double **) malloc (newSize * sizeof (double *));

    double **p1 = (double **) malloc (newSize * sizeof (double *));
    double **p2 = (double **) malloc (newSize * sizeof (double *));
    double **p3 = (double **) malloc (newSize * sizeof (double *));
    double **p4 = (double **) malloc (newSize * sizeof (double *));
    double **p5 = (double **) malloc (newSize * sizeof (double *));
    double **p6 = (double **) malloc (newSize * sizeof (double *));
    double **p7 = (double **) malloc (newSize * sizeof (double *));

    for (int i = 0; i < newSize; i++) 
    {
        a11[i] = (double *) malloc (newSize * sizeof (double));
        a12[i] = (double *) malloc (newSize * sizeof (double));
        a21[i] = (double *) malloc (newSize * sizeof (double));
        a22[i] = (double *) malloc (newSize * sizeof (double));

        b11[i] = (double *) malloc (newSize * sizeof (double));
        b12[i] = (double *) malloc (newSize * sizeof (double));
        b21[i] = (double *) malloc (newSize * sizeof (double));
        b22[i] = (double *) malloc (newSize * sizeof (double));

        s1[i] = (double *) malloc (newSize * sizeof (double));
        s2[i] = (double *) malloc (newSize * sizeof (double));
        s3[i] = (double *) malloc (newSize * sizeof (double));
        s4[i] = (double *) malloc (newSize * sizeof (double));

        t1[i] = (double *) malloc (newSize * sizeof (double));
        t2[i] = (double *) malloc (newSize * sizeof (double));
        t3[i] = (double *) malloc (newSize * sizeof (double));
        t4[i] = (double *) malloc (newSize * sizeof (double));

        p1[i] = (double *) malloc (newSize * sizeof (double));
        p2[i] = (double *) malloc (newSize * sizeof (double));
        p3[i] = (double *) malloc (newSize * sizeof (double));
        p4[i] = (double *) malloc (newSize * sizeof (double));
        p5[i] = (double *) malloc (newSize * sizeof (double));
        p6[i] = (double *) malloc (newSize * sizeof (double));
        p7[i] = (double *) malloc (newSize * sizeof (double));
    }
   for (int i = 0; i < newSize; i++) 
    {
        for (int j = 0; j < newSize; j++) 
        {
            a11[i][j] = A[i][j]; 
            a12[i][j] = A[i][j + newSize];
            a21[i][j] = A[i + newSize][j];
            a22[i][j] = A[i + newSize][j + newSize]; 

            b11[i][j] = B[i][j];
            b12[i][j] = B[i][j + newSize];
            b21[i][j] = B[i + newSize][j];
            b22[i][j] = B[i + newSize][j + newSize];
        }
    }

    s1 = matrix_add   (a21, a22, s1, newSize);
    s2 = matrix_subtr (s1, a11, s2, newSize);
    s3 = matrix_subtr (a11, a21, s3, newSize);
    s4 = matrix_subtr (a12, s2, s4, newSize);

    t1 = matrix_subtr (b12, b11, t1, newSize);
    t2 = matrix_subtr (b22, t1, t2, newSize);
    t3 = matrix_subtr (b22, b12, t3, newSize);
    t4 = matrix_subtr (b21, t2, t4, newSize);

    p1 = strassen(newSize, a11, b11, p1, c11, c12, c21, c22); 
    p2 = strassen(newSize, a12, b21, p2, c11, c12, c21, c22);
    p3 = strassen(newSize, s1, t1, p3  , c11, c12, c21, c22);
    p4 = strassen(newSize, s2, t2, p4  , c11, c12, c21, c22); 
    p5 = strassen(newSize, s3, t3, p5  , c11, c12, c21, c22);
    p6 = strassen(newSize, s4, b22, p6 , c11, c12, c21, c22);
    p7 = strassen(newSize, a22, t4, p7 , c11, c12, c21, c22);    

    // calculate c11, c12, c21, c22
    for (int i = 0; i < newSize; i++) 
    {
        for (int j = 0; j < newSize; j++) 
        {
            c11[i][j] = p1[i][j] + p2[i][j];
            c12[i][j] = p1[i][j] + p4[i][j] + p3[i][j] + p6[i][j]; 
            c21[i][j] = p1[i][j] + p4[i][j] + p5[i][j] + p7[i][j]; 
            c22[i][j] = p1[i][j] + p4[i][j] + p5[i][j] + p3[i][j];
        }
    }

    // join partitions back into C
    for (int i = 0; i < newSize; i++) 
    {
        for (int j = 0; j < newSize; j++) 
        {
            C[i][j] = c11[i][j];
            C[i][j + newSize] = c12[i][j];
            C[i + newSize][j] = c21[i][j];
            C[i + newSize][j + newSize] = c22[i][j]; 
        }
    }
    for (int i = 0; i < newSize; i++)
    {
        free (a11[i]); free (a12[i]); free (a21[i]); free (a22[i]);
        free (b11[i]); free (b12[i]); free (b21[i]); free (b22[i]);
        free (s1[i]); free (s2[i]); free (s3[i]); free (s4[i]);
        free (t4[i]); free (t3[i]); free (t2[i]); free (t1[i]);
    }
    // free memory
    free(a11); free(a12); free(a21); free(a22);
    free(b11); free(b12); free(b21); free(b22);
    free(s1); free(s2); free(s3); free(s4);
    free(t1); free(t2); free(t3); free(t4);

    return C;
}

struct pad_struct pad_matrices(double **a, double **b, int tile_size, int size) {
  struct pad_struct ret;
  int total_size = tile_size * 16;

  double **res_a = (double **) calloc(total_size, sizeof(double*));
  double **res_b = (double **) calloc(total_size, sizeof(double*));

  for (int i = 0; i < total_size; i++) {
    res_a[i] = (double *) calloc (total_size, sizeof(double));
    res_b[i] = (double *) calloc (total_size, sizeof(double));
  }

  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {

        res_a[i][j] = 1;
        res_b[i][j] = 2;
    }
  }

  ret.padded_a = res_a;
  ret.padded_b = res_b;
  ret.padded_size = total_size;
  ret.original_size = size;
  return ret;
}

double** unpad_matrix(double **m, int padded_size, int orig_size) {
    double **res = (double **) calloc(orig_size, sizeof(double*));
    for (int i = 0; i < orig_size; i++) {
        res[i] = (double *) calloc (orig_size, sizeof(double));
    }

    for (int i = 0; i < orig_size; i++) {
        for (int j = 0; j < orig_size; j++) {
            res[i][j] = m[i][j];
        }
    }
    return res;
}

bool power_of_two(int x) {
  return (x & (x - 1)) == 0;
}

void naive (double **A, double **B, double **res, int size)
{
    double r = 0.0;
    for (int i = 0; i < size; i++)
    {
        for (int k = 0; k < size; k++)
        {
            r = A[i][k];
            for (int j = 0; j < size; j++)
            {
                res[i][j] += r * B[k][j];
            }
        }
    }
}

double** matrix_add (double **A, double **B, double **res, int size)
{
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            res[i][j] = A[i][j] + B[i][j];
        }
    }

    return res;
}

double** matrix_subtr (double **A, double **B, double **res, int size)
{
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            res[i][j] = A[i][j] - B[i][j];
        }
    }

    return res;
}

void print_matrix (double **matrix, int size)
{
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            printf("%.1f ", matrix[i][j]);
        }
        printf("\n");
    }
}

void jki(double **a, double **b, double **c, int n) {
    double r = 0.0;
    for (int j = 0; j < n; j++) {
        for (int k = 0; k < n; k++) {
            r = b[k][j];
            for (int i = 0; i < n; i++) {
                c[i][j] += a[i][k] * r;
            }
        }
    }
}

void kji(double **a, double **b, double **c, int n) {
    double r = 0.0;
    for (int k = 0; k < n; k++) {
        for (int j = 0; j < n; j++) {
            r = b[k][j];
            for (int i = 0; i < n; i++) {
                c[i][j] += a[i][k] * r;
            }
        }
    }
}

void ijk(double **a, double **b, double **c, int n) {
    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            sum = 0.0;
            for (int k = 0; k < n; k++) {
                sum += a[i][k] * b[k][j];
            }
            c[i][j] += sum;
        }
    }
}

void jik(double **a, double **b, double **c, int n) {
    double sum = 0.0;
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < n; i++) {
            sum = 0.0;
            for (int k = 0; k < n; k++) {
                sum += a[i][k] * b[k][j];
            }
            c[i][j] += sum;
        }
    }
}

void kij(double **a, double **b, double **c, int n) {
    double r = 0.0;
    for (int k = 0; k < n; k++) {
        for (int i = 0; i < n; i++) {
            r = a[i][k];
            for (int j = 0; j < n; j++) {
                c[i][j] += r * b[k][j];
            }
        }
    }
}  

void ikj(double **a, double **b, double **c, int n) {
    double r = 0.0;
    for (int i = 0; i < n; i++) {
        for (int k = 0; k < n; k++) {
            r = a[i][k];
            for (int j = 0; j < n; j++) {
                c[i][j] += r * b[k][j];
            }
        }
    }
}