#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>
//#include <Accelerate/Accelerate.h>
//#include <omp.h>

struct pad_struct {
  double *padded_a;
  double *padded_b;
  int     padded_size;
  int     original_size;
};

double*  strassen  (int, double *, double *, double *, double *, 
              double *, double *, double *, double *, double *, 
              double *, double *, double *, double *, double *, 
              double *, double *, double *, double **);
void         morton_naive      (double *, double *, double *, uint32_t);
void         matrix_subtr      (double *, double *, double *, uint32_t);
void         matrix_add        (double *, double *, double *, uint32_t);
void         naive             (double *, double *, double *, uint32_t);
void         bijk(double *, double *, double *, int, int, int, int);
unsigned int S                 (unsigned int, unsigned int);
void         convertFromMorton (double *, double *, int);
void         convertToMorton   (double *, double *, int);
void         print_matrix      (double *, uint32_t);
void         spcl_print        (double *, int, int);
double       rand_from         (double, double);
int          L_r               (int, int, int);  
void         print_bits        (unsigned int);
int          layout            (int, int);
bool         power_of_two(int);
int          find_tilesize(int);
struct       pad_struct pad_matrices(double *, double *, int, int);
double*      unpad_matrix(double *, int, int);

const unsigned int E[] = {0x55555555, 0x33333333, 0x0F0F0F0F, 0x00FF00FF};
const unsigned int F[] = {1, 2, 4, 8};
unsigned int TR = 2, TC = 2;
int count = 0;
int final_idx = 0;

// int sizes[] = {32, 64, 128, 256, 512, 1024, 2048};
//int sizes[] = {100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200}; 
int sizes[] = {25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 375, 400, 425, 450, 475, 500, 525, 550, 575, 600, 625, 650, 675, 700};
int main() { 
  bool padded; 
  struct pad_struct padded_matrices; 
  for (int s = 0; s < 28; s++) {
        struct timespec tick, tock;
        long long elapsed_ns;
        // int n = sizes[s];
        int n = sizes[s];
        int orig_size = n;
        padded = false;
        double *A   = (double *) calloc (n * n, sizeof (double));
        double *B   = (double *) calloc (n * n, sizeof (double));
        double *E   = (double *) calloc (n * n, sizeof (double));
        for (int i = 0; i < n * n; i++) {
            A[i] = i;
            B[i] = i;
        }
        clock_gettime(CLOCK_REALTIME, &tick);
        int newSize;
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

        
        double *C   = (double *) calloc (n * n, sizeof (double));
        double *A_Morton =  (double *) calloc (n * n, sizeof (double));
        double *B_Morton =  (double *) calloc (n * n, sizeof (double));
        double *res = (double *) calloc (n * n, sizeof (double));
        double *D   = (double *) calloc (n * n, sizeof (double));

        double *p1  = (double *) calloc (newSize * newSize, sizeof (double));
        double *p2  = (double *) calloc (newSize * newSize, sizeof (double));
        double *p3  = (double *) calloc (newSize * newSize, sizeof (double));
        double *p4  = (double *) calloc (newSize * newSize, sizeof (double));
        double *p5  = (double *) calloc (newSize * newSize, sizeof (double));
        double *p6  = (double *) calloc (newSize * newSize, sizeof (double));
        double *p7  = (double *) calloc (newSize * newSize, sizeof (double));

        double *final_p1  = (double *) malloc (newSize * newSize * sizeof (double));
        double *final_p2  = (double *) malloc (newSize * newSize * sizeof (double));
        double *final_p3 = (double *) malloc (newSize * newSize * sizeof (double));
        double *final_p4 = (double *) malloc (newSize * newSize * sizeof (double));
        double *final_p5 = (double *) malloc (newSize * newSize * sizeof (double));
        double *final_p6 = (double *) malloc (newSize * newSize * sizeof (double));
        double *final_p7 = (double *) malloc (newSize * newSize * sizeof (double));
        double *final_p_addresses[7] = {final_p1, final_p2, final_p3, final_p4, final_p5, final_p6, final_p7};

        double *u1  = (double *) calloc (newSize * newSize, sizeof (double));
        double *u3  = (double *) calloc (newSize * newSize, sizeof (double));
        double *u2  = (double *) calloc (newSize * newSize, sizeof (double));
        double *u5  = (double *) calloc (newSize * newSize, sizeof (double));
        double *u4  = (double *) calloc (newSize * newSize, sizeof (double));
        double *u6  = (double *) calloc (newSize * newSize, sizeof (double));
        double *u7  = (double *) calloc (newSize * newSize, sizeof (double));

        // initialize A and B 
        if (padded) {
            convertToMorton(padded_matrices.padded_a, A_Morton, n);
            convertToMorton(padded_matrices.padded_b, B_Morton, n);
            strassen(n, A_Morton, B_Morton, C, p1, p2, p3, p4, p5, p6, p7, u1, u2, u3, u4, u5, u6, u7, final_p_addresses);
            C = unpad_matrix(C, padded_matrices.padded_size, padded_matrices.original_size);
        } else {
            convertToMorton(A, A_Morton, n);
            convertToMorton(B, B_Morton, n);
            strassen(n, A_Morton, B_Morton, C, p1, p2, p3, p4, p5, p6, p7, u1, u2, u3, u4, u5, u6, u7, final_p_addresses);
        }

        convertFromMorton(res, C, n);
       
        clock_gettime(CLOCK_REALTIME, &tock);
        elapsed_ns = (tock.tv_sec - tick.tv_sec) * 1000000000LL + (tock.tv_nsec - tick.tv_nsec);
        printf("For size %d, Strassen CLOCK_REALTIME Clock elapsed: %lld ns\n", orig_size, elapsed_ns);
        A = (double *) calloc (n * n, sizeof (double));
        B = (double *) calloc (n * n, sizeof (double));
        for (int i = 0; i < n * n; i++) {
            A[i] = i;
            B[i] = i;
        }
        clock_gettime(CLOCK_REALTIME, &tick);
        morton_naive(A, B, D, n);
        clock_gettime(CLOCK_REALTIME, &tock);
        elapsed_ns = (tock.tv_sec - tick.tv_sec) * 1000000000LL + (tock.tv_nsec - tick.tv_nsec);
        printf("For size %d, Naive CLOCK_REALTIME elapsed: %lld ns\n", orig_size, elapsed_ns);
         
        //clock_gettime(CLOCK_REALTIME, &tick);
        //bijk(A, B, D, orig_size, orig_size, orig_size, 10);
        //clock_gettime(CLOCK_REALTIME, &tock);
        //elapsed_ns = (tock.tv_sec - tick.tv_sec) * 1000000000LL + (tock.tv_nsec - tick.tv_nsec);
        //printf("For size %d, BIJK CLOCK_REALTIME elapsed: %lld ns\n", orig_size, elapsed_ns);

        free(A); free(B); free(C); free(A_Morton); free(B_Morton); free(res); 
        free(p1); free(p2); free(p3); free(p4); free(p5); free(p6); free(p7);
        free(u1); free(u2); free(u3); free(u4); free(u5); free(u6); free(u7);
    }
    return 0;
}

double* strassen (int size, double *A, double *B, double *C,
                  double *p1, double *p2, double *p3, 
                  double *p4, double *p5, double *p6, double *p7, double *u1, 
                  double *u2, double *u3, double *u4, double *u5, double *u6, 
                  double *u7, double **final_p_addresses) {
    
    if (size == TR) {
        // double r = 0.0;
        // for (uint32_t i = 0; i < size; i++) {
        //     for (uint32_t k = 0; k < size; k++) {
        //         r = A[i * size + k];
        //         printf("A[%d]=%.1f\n", i * size + k, A[i * size + k]);
        //         for (uint32_t j = 0; j < size; j++) {
        //             C[i * size + j] += r * B[k * size + j];
        //             printf("B[%d]=%.1f\n", k * size + j, B[k * size + j]);
        //             printf("C[%d] = %.1f\n", i * size + j, C[i * size + j]);
        //         }
        //     }
        // }
        // printf("multiplying: \n");
        // print_matrix(A, size);
        // printf("with:\n");
        // print_matrix(B, size);
        // printf("in multiplying\n");
        morton_naive(A, B, C, size);
        // printf("multiplied\n");
        // print_matrix(A, size);
        // printf("by\n");
        // print_matrix(B, size);
        // printf("result is \n");
        // print_matrix(C, size);
//        C[0] = A[0] * B[0] + A[1] * B[2];
//        C[1] = A[0] * B[1] + A[1] * B[3];
//       C[2] = A[2] * B[0] + A[3] * B[2];
//        C[3] = A[2] * B[1] + A[3] * B[3];
        // printf("multiplying %f with %f resulting in %f\n", A[0], B[0], C[0]);
        return C;
    }
    
    int newSize = size / 2;
    // printf("__A is__ \n");
    // print_matrix(A, size);
    // printf("__B__ is \n");
    // print_matrix(B, size);
     
    // // A[0:15]
    // printf("**a11** is \n");
    // spcl_print(A, 0, newSize);
    // // A[16:31]
    // printf("\n**a12** is \n");
    // spcl_print(A, (newSize * newSize), newSize);
    
    // // A[32:47]
    // printf("\n**a21** is \n");
    // spcl_print(A, (newSize * newSize) * 2, newSize);
    
    // //A[48:64]
    // printf("\n**a22** is \n");
    // spcl_print(A, (newSize * newSize) * 3, newSize);

    // printf("\n**b11** is \n");
    // spcl_print(B, 0, newSize);
    // // A[16:31]
    // printf("\n**b12** is \n");
    // spcl_print(B, (newSize * newSize), newSize);
    
    // // A[32:47]
    // printf("\n**b21** is \n");
    // spcl_print(B, (newSize * newSize) * 2, newSize);
    
    // //A[48:64]
    // printf("\n**b22** is \n");
    // spcl_print(B, (newSize * newSize) * 3, newSize);
    // a11 = A + 0
    // a12 = A + ((newSize * newSize) * 1)
    // a21 = A + ((newSize * newSize) * 2)
    // a22 = A + ((newSize * newSize) * 3)
    double *S  = (double *) calloc (size * size, sizeof (double));
    double *T  = (double *) calloc (size * size, sizeof (double));

    // S1 = A21 + A22
    matrix_add   (A + ((newSize * newSize) * 2) , A + ((newSize * newSize) * 3), S, newSize);
    // // S2 = S1 - A11
    matrix_subtr (S, A,  S + ((newSize * newSize) * 1), newSize);
    // // S3 = A11 - A21
    matrix_subtr (A, A + ((newSize * newSize) * 2), S + ((newSize * newSize) * 2), newSize);
    // // S4 = A12 - S2
    matrix_subtr (A + ((newSize * newSize) * 1), S + ((newSize * newSize) * 1),  S + ((newSize * newSize) * 3), newSize);

    // printf("\n**s1** is \n");
    // print_matrix(S, newSize);
    // printf("**s2** is \n");
    // print_matrix(S + ((newSize * newSize) * 1), newSize);
    // printf("**s3** is \n");
    // print_matrix(S + ((newSize * newSize) * 2), newSize);
    // printf("**s4** is \n");
    // print_matrix(S + ((newSize * newSize) * 3), newSize);

    // T1 = B12 - B11
    matrix_subtr (B + ((newSize * newSize) * 1), B, T, newSize);
    // // T2 = B22 - T1
    matrix_subtr (B + ((newSize * newSize) * 3), T,  T + ((newSize * newSize) * 1), newSize);
    // // T3 = B22 - B12
    matrix_subtr (B + ((newSize * newSize) * 3), B + ((newSize * newSize) * 1), T + ((newSize * newSize) * 2), newSize);
    // // T4 = B21 - T2
    matrix_subtr (B + ((newSize * newSize) * 2), T + ((newSize * newSize) * 1),  T + ((newSize * newSize) * 3), newSize);

    // printf("**t1** is \n");
    // print_matrix(T, newSize);
    // printf("**t2** is \n");
    // print_matrix(T + ((newSize * newSize) * 1), newSize);
    // printf("**t3** is \n");
    // print_matrix(T + ((newSize * newSize) * 2), newSize);
    // printf("**t4** is \n");
    // print_matrix(T + ((newSize * newSize) * 3), newSize);

    // P1 = A11 * B11
    // printf("********** calling strassen for P1 **********\n");
    p1 = strassen(newSize, A, B, p1, p1, p2, p3, p4, p5, p6, p7, u1, u2, u3, u4, u5, u6, u7, final_p_addresses); 
    
    // P2 = A12 * B21
    // printf("********** calling strassen for P2 **********\n");
    p2 = strassen(newSize, A + ((newSize * newSize) * 1), B + ((newSize * newSize) * 2), p2, p1, p2, p3, p4, p5, p6, p7, 
    u1, u2, u3, u4, u5, u6, u7, final_p_addresses);

    // P3 = S1 * T1
    // printf("********** calling strassen for P3 **********\n");
    p3 = strassen(newSize, S, T, p3, p1, p2, p3, p4, p5, p6, p7, u1, u2, u3, u4, u5, u6, u7, final_p_addresses);

    // P4 = S2 * T2
    // printf("********** calling strassen for P4 **********\n");
    p4 = strassen(newSize, S + ((newSize * newSize) * 1), T + ((newSize * newSize) * 1), p4, p1, p2, p3, p4, p5, p6, p7, u1, u2, u3, u4, u5, u6, u7, final_p_addresses); 
    
    // P5 = S3 * T3
    // printf("********** calling strassen for P5 **********\n");
    p5 = strassen(newSize, S + ((newSize * newSize) * 2), T + ((newSize * newSize) * 2), p5, p1, p2, p3, p4, p5, p6, p7, u1, u2, u3, u4, u5, u6, u7, final_p_addresses);
    
    // P6 = S4 * B22
    // printf("********** calling strassen for P6 **********\n");
    p6 = strassen(newSize, S + ((newSize * newSize) * 3), B + ((newSize * newSize) * 3), p6, p1, p2, p3, p4, p5, p6, p7, u1, u2, u3, u4, u5, u6, u7, final_p_addresses);
    
    // P7 = A22 * T4
    // printf("********** calling strassen for P7 **********\n");
    p7 = strassen(newSize, A + ((newSize * newSize) * 3), T + ((newSize * newSize) * 3), p7, p1, p2, p3, p4, p5, p6, p7, u1, u2, u3, u4, u5, u6, u7, final_p_addresses);    

    // printf("**p1** is \n");
    // print_matrix(p1, newSize);
    // printf("**p2** is \n");
    // print_matrix(p2, newSize);
    // printf("**p3** is \n");
    // print_matrix(p3, newSize);
    // printf("**p4** is \n");
    // print_matrix(p4, newSize);
    // printf("**p5** is \n");
    // print_matrix(p5, newSize);
    // printf("**p6** is \n");
    // print_matrix(p6, newSize);
    // printf("**p7** is \n");
    // print_matrix(p7, newSize);
    // printf("count is %d\n", count);

    if (count < 7) {
        for (int i = 0 ; i < newSize * newSize; i++) {
            u1[i] = p1[i] + p2[i];
            u2[i] = p1[i] + p4[i];
            u3[i] = u2[i] + p5[i];
            u4[i] = u3[i] + p7[i];
            u5[i] = u3[i] + p3[i];
            u6[i] = u2[i] + p3[i];
            u7[i] = u6[i] + p6[i];
        }
    } else {
        // printf("in else!");
        double *finalp1 = final_p_addresses[0];
        double *finalp2 = final_p_addresses[1];
        double *finalp3 = final_p_addresses[2];
        double *finalp4 = final_p_addresses[3];
        double *finalp5 = final_p_addresses[4];
        double *finalp6 = final_p_addresses[5];
        double *finalp7 = final_p_addresses[6];
        for (int i = 0 ; i < newSize * newSize; i++) {
            u1[i] = finalp1[i] + finalp2[i];
            u2[i] = finalp1[i] + finalp4[i];
            u3[i] = u2[i] + finalp5[i];
            u4[i] = u3[i] + finalp7[i];
            u5[i] = u3[i] + finalp3[i];
            u6[i] = u2[i] + finalp3[i];
            u7[i] = u6[i] + finalp6[i];
        }
        count = 0;
        final_idx = 0;
    }
    
    // join partitions back into C
    // check these conditions
    for (int i = 0; i < newSize * newSize; i++) {
        C[i + 0] = u1[i];
        C[i + ((newSize * newSize) * 1)] = u4[i];
        C[i + ((newSize * newSize) * 2)] = u7[i];
        C[i + ((newSize * newSize) * 3)] = u5[i];
    }

    if (newSize == TR) { //if newSize == TR
        // printf("final!. size is %d\n", size);
        double* cur_p = final_p_addresses[final_idx];
        for (int i = 0; i < size * size; i++) {
            // printf("look\n");
            cur_p[i] = C[i];
            // printf("after\n");
        }
        // printf("final_p_addresses[%d] is:\n", final_idx);
        // print_matrix(final_p_addresses[final_idx], size);
        final_idx++;
        count++;
        // printf("count is %d\n", count);
    }
    free(S);
    free(T);
    // printf("at the end, C is \n");  
    // print_matrix(C, size);
    return C;
}

double *pad_matrix(double *a, int tile_size, int size) {
  int total_size = tile_size + size;
  double *res = (double *) calloc(total_size * total_size, sizeof(double));
  for (int i = 0; i < size; i++) {
    res[i] = a[i];
  }

  return  res;
}

// helper method to pad matrices that are not of size powers of 2.
struct pad_struct pad_matrices(double *a, double *b, int tile_size, int size) {
  struct pad_struct ret;
  int total_size = tile_size * 16;
  double *res_a = (double *) calloc(total_size * total_size, sizeof(double));
  double *res_b = (double *) calloc(total_size * total_size, sizeof(double));

  for (int i = 0; i < size * size; i++) {
    res_a[i] = a[i];
    res_b[i] = b[i];
  }
  
  ret.padded_a = res_a;
  ret.padded_b = res_b;
  ret.padded_size = total_size;
  ret.original_size = size;
  return ret;
}

double* unpad_matrix(double *m, int padded_size, int orig_size) {
    double *res = (double *) calloc(orig_size * orig_size, sizeof(double));
    for (int i = 0; i < orig_size * orig_size; i++) {
        res[i] = m[i];
    }
    
    return res;
}

// helper method to check if a number is a power of two.
bool power_of_two(int x) {
  return (x & (x - 1)) == 0;
}

int find_tilesize(int size) {
  // example: square matrix size of 513. Chosing size of 32 would mess it up, but choosing size of 33 
  // would make it more efficient.
  // Would add 15 of padding, making size = 528. 528 / 2 / 2 / 2 / 2 == 33. 

  int initial_size = 32;
  int tile_size = initial_size;

  while ((tile_size *= 2) < size);
  while (tile_size > size * 1.1) {
    initial_size += 1;
    tile_size = initial_size;
    while ((tile_size *= 2) < size);
  }

  return initial_size;
}

// source: https://dl.acm.org/doi/pdf/10.1145/305619.305645 pg. 225
int layout (int i, int j)
{
    // i == row, j == col. 
    // ti = T (i, tr) = i / tr
    // tj = T (j, tc) = j / tc
    // fi = F (i, tr) = i % tr
    // fj = F (j, tc) = j % tc
    // TR * TC * S (ti, tj) + Lr (fi, fj, tr, tc)
    return TR * TC * S (i / TR, j / TC) + L_r (i % TR, j % TC, TC);

}

int L_r (int i, int j, int n)
{
    return n * i + j;
}

void convertToMorton(double *matrix, double *morton, int size)
{
    for (int row = 0; row < size; row++)
    {
        for (int col = 0; col < size; col++)
        {
            int res = layout (row, col);
            morton[res] = matrix[row * size + col];
        }
    }
}

void convertFromMorton(double *matrix, double *morton, int size) {
    for (int row = 0; row < size; row++) {
        for (int col = 0; col < size; col++) {
            int res = layout(row, col);
            matrix[res] = morton[row * size + col];
        }
    }
}

// source: https://graphics.stanford.edu/~seander/bithacks.html#InterleaveBMN
unsigned int S (unsigned int x, unsigned int y)
{
    unsigned int z = 0;
    x = (x | (x << F[3])) & E[3];
    x = (x | (x << F[2])) & E[2];
    x = (x | (x << F[1])) & E[1];
    x = (x | (x << F[0])) & E[0];

    y = (y | (y << F[3])) & E[3];
    y = (y | (y << F[2])) & E[2];
    y = (y | (y << F[1])) & E[1];
    y = (y | (y << F[0])) & E[0];

    z = y | (x << 1);
    return z;
}

void matrix_subtr (double *a, double *b, double *res, uint32_t size)
{
#pragma omp parallel for 
    for (uint32_t i = 0; i < size * size; i++) {
        res[i] = a[i] - b[i];
    }
}

void matrix_add (double *a, double *b, double *res, uint32_t size)
{
#pragma omp parallel for
    for (uint32_t i = 0; i < size * size; i++) {
        res[i] = a[i] + b[i];
    }
}

void naive (double *a, double *b, double *res, uint32_t size) {
    double sum = 0.0;
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            sum = 0.0;
            for (int k = 0; k < size; k++) {
                sum += a[i * size + k] * b[k * size + j];
            }
            res[i * size + j] += sum;
        }
    }
}

void bijk(double *A, double *B, double *C, int m, int n, int p, int block) {
  for (int i = 0; i < m; i += block) {
    for (int k = 0; k < n; k += block) {
      for (int j = 0; j < p; j += block) {
        for (size_t b_i = i; b_i < i+block && b_i < m; b_i++) {
          for (size_t b_k = k; b_k < k+block && b_k < n; b_k++) {
            for (size_t b_j = j; b_j < j+block && b_j < p; b_j++) { 
              C[b_i*n + b_j] += A[b_i*n + b_k] * B[b_k*p + b_j];
            }
          }
        }
      }
    }
  }
}

void morton_naive (double *a, double *b, double *res, uint32_t size) {
    double r = 0.0;
// #pragma omp parallel for
    for (uint32_t i = 0; i < size; i++) {
        for (uint32_t k = 0; k < size; k++) {
            r = a[i * size + k];
            for (uint32_t j = 0; j < size; j++) {
              res[i * size + j] += r * b[k * size + j];
            }
        }
    }
}

double rand_from (double min, double max)
{
    double rand_decimal = ((double) rand() / RAND_MAX) * (max - min) + min;
    
    return (uint32_t) (rand_decimal * 100.0 + 0.5) / 100.0;
}

void spcl_print (double *matrix, int starting, int size) {
    for (int i = 0; i < size * size; i++) {
        if (i % size == 0) {
            printf("\n");
        }
        printf("%.1f ", matrix[i + starting]);
    }
}

void print_matrix (double *matrix, uint32_t size)
{
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            printf("%.1f ", matrix[i * size + j]);
        }
        printf("\n");
    }
}

void print_bits (unsigned int num)
{
    unsigned int size = sizeof (unsigned int);
    unsigned int factor = 1 << (8 * size - 1);
    int iter = 0;
    while (factor)
    {
        printf("%i", factor & num ? 1 : 0);
        factor >>= 1;
        iter++;
        if (iter % 8 == 0)
        {
            printf(" ");
        }
    }
    printf("\n");
}
