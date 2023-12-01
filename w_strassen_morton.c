#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>
#include <papi.h>

struct pad_struct {
  double *padded_a;
  double *padded_b;
} pad_struct;

void  strassen  (int, double *, double *, double *, double *, 
              double *, double *, double *, double *, double *, 
              double *, double *, double *, double *, double *, 
              double *, double *, double *);
void         morton_naive      (double *, double *, double *, uint32_t);
void         matrix_subtr      (double *, double *, double *, uint32_t);
void         matrix_add        (double *, double *, double *, uint32_t);
void         naive             (double *, double *, double *, uint32_t);
unsigned int S                 (unsigned int, unsigned int);
void         convertFromMorton (double *, double *, int);
void         convertToMorton   (double *, double *, double *, double *, int);
void         print_matrix      (double *, uint32_t);
void         spcl_print        (double *, int, int);
double       rand_from         (double, double);
int          L_r               (int, int, int);  
void         print_bits        (unsigned int);
int          layout            (int, int);
bool power_of_two(int);
int find_tilesize(int);
// struct pad_struct *pad_matrices(double *, double *, int, int);
double *pad_matrix(double *, int, int);
double *unpad_matrix(double *, int, int);
void print_data(long long *);

const unsigned int E[] = {0x55555555, 0x33333333, 0x0F0F0F0F, 0x00FF00FF};
const unsigned int F[] = {1, 2, 4, 8};
unsigned int TR = 16, TC = 16;
int PAPI_one[4] = {
    PAPI_L1_DCA,
    PAPI_L1_DCM,
    PAPI_L2_DCA,
    PAPI_L2_DCM,
};

int sizes[] = {2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048};
int main() {
    int retval;
    if ((retval = PAPI_library_init(PAPI_VER_CURRENT)) == -1) {
        fprintf(stderr, "Error initializing PAPI! %s\n", PAPI_strerror(retval));
        return -1;
    }
    long long counters[4];
    for (int s = 0; s < 11; s++) {
        struct timespec tick, tock;
        long long elapsed_ns;
        int n = sizes[s];
        int newSize = n / 2;
        double *A   = (double *) calloc (n * n, sizeof (double));
        double *B   = (double *) calloc (n * n, sizeof (double));
        double *C   = (double *) calloc (n * n, sizeof (double));
        // initialize A and B 
        for (int i = 0; i < n * n; i++) {
            A[i] = i;
            B[i] = i;
        }
        // start timing from here.
        // if (PAPI_start_counters(PAPI_one, 4) != 0) {
        //     fprintf(stderr, "An error has occurred\n");
        // }
        double *A_Morton =  (double *) calloc (n * n, sizeof (double));
        double *B_Morton =  (double *) calloc (n * n, sizeof (double));
        double *res = (double *) calloc (n * n, sizeof (double));
        double *p1  = (double *) calloc (newSize * newSize, sizeof (double));
        double *p3  = (double *) calloc (newSize * newSize, sizeof (double));
        double *p2  = (double *) calloc (newSize * newSize, sizeof (double));
        double *p5  = (double *) calloc (newSize * newSize, sizeof (double));
        double *p4  = (double *) calloc (newSize * newSize, sizeof (double));
        double *p6  = (double *) calloc (newSize * newSize, sizeof (double));
        double *p7  = (double *) calloc (newSize * newSize, sizeof (double));
        double *u1  = (double *) calloc (newSize * newSize, sizeof (double));
        double *u3  = (double *) calloc (newSize * newSize, sizeof (double));
        double *u2  = (double *) calloc (newSize * newSize, sizeof (double));
        double *u5  = (double *) calloc (newSize * newSize, sizeof (double));
        double *u4  = (double *) calloc (newSize * newSize, sizeof (double));
        double *u6  = (double *) calloc (newSize * newSize, sizeof (double));
        double *u7  = (double *) calloc (newSize * newSize, sizeof (double));
        clock_gettime(CLOCK_REALTIME, &tick);
        
        // if (n % 4 != 0) {
        //     int tile_size = find_tilesize(n);
        //     double *padded_a = pad_matrix(A, tile_size, n);
        //     double *padded_b = pad_matrix(B, tile_size, n);
        //     convertToMorton(padded_a, A_Morton, n + tile_size);
        //     convertToMorton(padded_b, B_Morton, n + tile_size);
        //     double *C = (double *) calloc ((n + tile_size) * (n + tile_size), sizeof(double));
        //     strassen(n, A_Morton, B_Morton, C, p1, p2, p3, p4, p5, p6, p7, u1, u2, u3, u4, u5, u6, u7);
        //     // C = unpad_matrix(C, tile_size, n);
        // } else {
        //     if (TR > n) {
        //         naive(A, B, C, n);
        //     } else {
        //         convertToMorton(A, A_Morton, n);
        //         convertToMorton(B, B_Morton, n);
        //         double *C = (double *) calloc (n * n, sizeof(double));
        //         strassen(n, A_Morton, B_Morton, C, p1, p2, p3, p4, p5, p6, p7, u1, u2, u3, u4, u5, u6, u7);
        //         convertFromMorton(res, C, n);
        //     }
        // }
        if (TR > n) {
            naive(A, B, C, n);
        } else {
            convertToMorton(A, B, A_Morton, B_Morton, n);
            
            strassen(n, A_Morton, B_Morton, C, p1, p2, p3, p4, p5, p6, p7, u1, u2, u3, u4, u5, u6, u7);
            convertFromMorton(res, C, n);
        }
        
        // if (PAPI_stop_counters(counters, 4) != 0) {
        //     fprintf(stderr, "An error has occured\n");
        //     return -1;
        // }
        // printf("\033[1;32mTiming data for size=%d\033[0m\n", n);
        // print_data(counters);

        clock_gettime(CLOCK_REALTIME, &tock);
        elapsed_ns = (tock.tv_sec - tick.tv_sec) * 1000000000LL + (tock.tv_nsec - tick.tv_nsec);

        printf("For size %d, CLOCK_REALTIME Clock elapsed: %lld ns\n", n, elapsed_ns);
        free(A); free(B); free(C); free(A_Morton); free(B_Morton); free(res); 
        free(p1); free(p2); free(p3); free(p4); free(p5); free(p6); free(p7);
        free(u1); free(u2); free(u3); free(u4); free(u5); free(u6); free(u7);
    }
    return 0;
}

void strassen (int size, double *A, double *B, double *C,
                  double *p1, double *p2, double *p3, 
                  double *p4, double *p5, double *p6, double *p7, double *u1, 
                  double *u2, double *u3, double *u4, double *u5, double *u6, 
                  double *u7) {
    if (size == TR) {
        morton_naive(A, B, C, size);
        return;
    }
    
    int newSize = size / 2;
    // printf("__A is__ \n");
    // print_matrix(A, size);
    // printf("__B__ is \n");
    // print_matrix(B, size);
     
    // A[0:15]
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

    // printf("**b11** is \n");
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

    // double *s1 = (double *) calloc (newSize * newSize, sizeof (double));
    // double *s2 = (double *) calloc (newSize * newSize, sizeof (double));
    // double *s3 = (double *) calloc (newSize * newSize, sizeof (double));
    // double *s4 = (double *) calloc (newSize * newSize, sizeof (double));

    // double *t1 = (double *) calloc (newSize * newSize, sizeof (double));
    // double *t2 = (double *) calloc (newSize * newSize, sizeof (double));
    // double *t3 = (double *) calloc (newSize * newSize, sizeof (double));
    // double *t4 = (double *) calloc (newSize * newSize, sizeof (double));

    // matrix_add   (A + ((newSize * newSize) * 2) , A + ((newSize * newSize) * 3), s1, newSize);
    // // // S2 = S1 - A11
    // matrix_subtr (s1, A,  s2, newSize);
    // // // S3 = A11 - A21
    // matrix_subtr (A, A + ((newSize * newSize) * 2), s3, newSize);
    // // // S4 = A12 - S2
    // matrix_subtr (A + ((newSize * newSize) * 1), s2,  s4, newSize);

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

    //     // T1 = B12 - B11
    // matrix_subtr (B + ((newSize * newSize) * 1), B, t1, newSize);
    // // // T2 = B22 - T1
    // matrix_subtr (B + ((newSize * newSize) * 3), t1,  t2, newSize);
    // // // T3 = B22 - B12
    // matrix_subtr (B + ((newSize * newSize) * 3), B + ((newSize * newSize) * 1), t3, newSize);
    // // // T4 = B21 - T2
    // matrix_subtr (B + ((newSize * newSize) * 2), S + ((newSize * newSize) * 1),  t4, newSize);

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
    strassen(newSize, A, B, p1, p1, p2, p3, p4, p5, p6, p7, u1, u2, u3, u4, u5, u6, u7); 
    // P2 = A12 * B21
    // printf("********** calling strassen for P2 **********\n");
    strassen(newSize, A + ((newSize * newSize) * 1), B + ((newSize * newSize) * 2), p2, p1, p2, p3, p4, p5, p6, p7, 
    u1, u2, u3, u4, u5, u6, u7);
    // P3 = S1 * T1
    // printf("********** calling strassen for P3 **********\n");
    strassen(newSize, S, T, p3, p1, p2, p3, p4, p5, p6, p7, u1, u2, u3, u4, u5, u6, u7);
    // P4 = S2 * T2
    // printf("********** calling strassen for P4 **********\n");
    strassen(newSize, S + ((newSize * newSize) * 1), T + ((newSize * newSize) * 1), p4, p1, p2, p3, p4, p5, p6, p7, u1, u2, u3, u4, u5, u6, u7); 
    // P5 = S3 * T3
    // printf("********** calling strassen for P5 **********\n");
    strassen(newSize, S + ((newSize * newSize) * 2), T + ((newSize * newSize) * 2), p5, p1, p2, p3, p4, p5, p6, p7, u1, u2, u3, u4, u5, u6, u7);
    // P6 = S4 * B22
    // printf("********** calling strassen for P6 **********\n");
    strassen(newSize, S + ((newSize * newSize) * 3), B + ((newSize * newSize) * 3), p6, p1, p2, p3, p4, p5, p6, p7, u1, u2, u3, u4, u5, u6, u7);
    // P7 = A22 * T4
    // printf("********** calling strassen for P7 **********\n");
    strassen(newSize, A + ((newSize * newSize) * 3), T + ((newSize * newSize) * 3), p7, p1, p2, p3, p4, p5, p6, p7, u1, u2, u3, u4, u5, u6, u7);    

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
    
    for (int i = 0 ; i < newSize * newSize; i++) {
        u1[i] = p1[i] + p2[i];
        u2[i] = p1[i] + p4[i];
        u3[i] = u2[i] + p5[i];
        u4[i] = u3[i] + p7[i];
        u5[i] = u3[i] + p3[i];
        u6[i] = u2[i] + p3[i];
        u7[i] = u6[i] + p6[i];
    }
    // join partitions back into C
    // check these conditions
    for (int i = 0; i < newSize * newSize; i++) {
        C[i + 0] = u1[i];
        C[i + ((newSize * newSize) * 1)] = u4[i];
        C[i + ((newSize * newSize) * 2)] = u7[i];
        C[i + ((newSize * newSize) * 3)] = u5[i];
    }
    
    free(S);
    free(T);
    // printf("at the end, C is \n");
    // print_matrix(C, size);
}

// helper method to pad matrix given the total tile size.
double *pad_matrix(double *source, int orig_size, int padded_size) {
    double *padded = (double *) calloc(padded_size * padded_size, sizeof(double));
    int padded_idx = 0;
    int diff = padded_size - orig_size;
    for (int i = 0; i < orig_size * orig_size; i++) {
        if (i > 0 && i % orig_size == 0) {
            padded_idx += diff;
        }
        padded[padded_idx] = source[i];
        padded_idx++;
    }
    
    return padded;
}

// helper method to pad matrices that are not of size powers of 2.
struct pad_struct *pad_matrices(double *a, double *b, int tile_size, int size) {
  struct pad_struct *ret;
  int total_size = tile_size + size;
  double *res_a = calloc(total_size * total_size, sizeof(double));
  double *res_b = calloc(total_size * total_size, sizeof(double));
  for (int i = 0; i < size; i++) {
    res_a[i] = a[i];
    res_b[i] = b[i];
  }
  ret->padded_a = res_a;
  ret->padded_b = res_b;

  return ret;
}

// helper method to check if a number is a power of two.
bool power_of_two(int x) {
  return (x & (x - 1)) == 0;
}

int find_tilesize(int size) {
  // example: square matrix size of 513. Chosing size of 32 would mess it up, but choosing size of 33 
  // would make it more efficient.
  // Would add 15 of padding, making size = 528. 528 / 2 / 2 / 2 / 2 == 33. 

  int initial_size = TR;
  int tile_size = initial_size;

  while ((tile_size *= 2) < size);
  while (tile_size > size * 1.1 && tile_size < 64) {
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

void convertToMorton(double *a, double *b, double *a_morton, double *b_morton, int size)
{
    for (int row = 0; row < size; row++)
    {
        for (int col = 0; col < size; col++)
        {
            int res = layout (row, col);
            a_morton[res] = a[row * size + col];
            b_morton[res] = b[row * size + col];
        }
    }
}

void print_data(long long *counters) {
    printf("L1 Miss Rate is: %0.3f\n", ((float) counters[1] / (float) counters[0]));
    printf("L2 Miss Rate is: %0.3f\n", ((float) counters[3] / (float) counters[2]));
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
    for (uint32_t i = 0; i < size * size; i++) {
        res[i] = a[i] - b[i];
    }
    // for (uint32_t i = 0; i < size * size; i++)
    // {
    //     for (uint32_t j = 0; j < size; j++)
    //     {
    //         res[i * size + j] = a[i * size + j] - b[i * size + j];
    //     }
    // }
}

void matrix_add (double *a, double *b, double *res, uint32_t size)
{
    for (uint32_t i = 0; i < size * size; i++) {
        res[i] = a[i] + b[i];
    }
    // for (uint32_t i = 0; i < size; i++)
    // {
    //     for (uint32_t j = 0; j < size; j++)
    //     {
    //         res[i * size + j] = a[i * size + j] + b[i * size + j];
    //     }
    // }
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

void morton_naive (double *a, double *b, double *res, uint32_t size) {
    double sum = 0.0;
    for (uint32_t i = 0; i < size; i++)
    {
        for (uint32_t j = 0; j < size; j++)
        {
            sum = 0.0;
            for (uint32_t k = 0; k < size; k++)
            {
              sum += a[i * size + k] * b[k * size + j]; 
            }
            res[i * size + j] = sum;
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
