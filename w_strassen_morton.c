#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>

void  strassen(int, double *, double *, double *, double *, 
              double *, double *, double *, double *, double *,
              double *, double *, double *, double *, double *,
              double *, double *, double *, double *, double *, 
              double *, double *, double *, double *, double *, 
              double *p6, double *p7);
void         matrix_subtr    (double *, double *, double *, uint32_t);
void         matrix_add      (double *, double *, double *, uint32_t);
void         naive           (double *, double *, double *, uint32_t);
unsigned int S               (unsigned int, unsigned int);
void         convertToMorton (double *, double *, int);
void         print_matrix    (double *, uint32_t);
double       rand_from       (double, double);
int          L_r             (int, int, int);  
void         print_bits      (unsigned int);
int          layout          (int, int);

const unsigned int E[] = {0x55555555, 0x33333333, 0x0F0F0F0F, 0x00FF00FF};
const unsigned int F[] = {1, 2, 4, 8};
const unsigned int TR = 2, TC = 2;

int main()
{
    int sizes[] = {8};
    for (int s = 0; s < 1; s++)
    {
        int n = sizes[s];
        int newSize = n / 2;
        // index into A as A[i * cur_size + j]
        double *A   = (double *) malloc (n * n * sizeof (double));
        double *B   = (double *) malloc (n * n * sizeof (double));
        double *A_Morton =  (double *) malloc (n * n * sizeof (double));
        double *B_Morton =  (double *) malloc (n * n * sizeof (double));
        double *C   = (double *) malloc (n * n * sizeof (double));
        double *D   = (double *) malloc (n * n * sizeof (double));
        
        double *a11 = (double *) malloc (newSize * newSize * sizeof (double));
        double *a12 = (double *) malloc (newSize * newSize * sizeof (double)); 
        double *a21 = (double *) malloc (newSize * newSize * sizeof (double));
        double *a22 = (double *) malloc (newSize * newSize * sizeof (double));

        double *b11 = (double *) malloc (newSize * newSize * sizeof (double));
        double *b12 = (double *) malloc (newSize * newSize * sizeof (double)); 
        double *b21 = (double *) malloc (newSize * newSize * sizeof (double));
        double *b22 = (double *) malloc (newSize * newSize * sizeof (double));

        double *s1 = (double *) malloc (newSize * newSize * sizeof (double));
        double *s2 = (double *) malloc (newSize * newSize * sizeof (double));
        double *s3 = (double *) malloc (newSize * newSize * sizeof (double));
        double *s4 = (double *) malloc (newSize * newSize * sizeof (double));

        double *t1 = (double *) malloc (newSize * newSize * sizeof (double));
        double *t2 = (double *) malloc (newSize * newSize * sizeof (double));
        double *t3 = (double *) malloc (newSize * newSize * sizeof (double));
        double *t4 = (double *) malloc (newSize * newSize * sizeof (double));

        double *p1 = (double *) malloc (newSize * newSize * sizeof (double));
        double *p3 = (double *) malloc (newSize * newSize * sizeof (double));
        double *p2 = (double *) malloc (newSize * newSize * sizeof (double));
        double *p5 = (double *) malloc (newSize * newSize * sizeof (double));
        double *p4 = (double *) malloc (newSize * newSize * sizeof (double));
        double *p7 = (double *) malloc (newSize * newSize * sizeof (double));
        double *p6 = (double *) malloc (newSize * newSize * sizeof (double));


        // initialize A and B 
        for (int i = 0; i < n * n; i++) 
        {
            A[i] = i;
            B[i] = i;
            // for (int j = 0; j < n; j++) 
            // {
            //     A[i * n + j] = i + j;
            //     B[i * n + j] = j;
            // }
        }

        convertToMorton(A, A_Morton, n);
        printf("A is: \n");
        print_matrix (A, n);
        printf("Initially A_Morton is:\n");
        print_matrix (A_Morton, n);
        convertToMorton(B, B_Morton, n);
        printf("B is: \n");
        print_matrix (B, n);
        printf("Initially B_Morton is:\n");
        print_matrix (B_Morton, n);
        // strassen (n, A_Morton, B_Morton, C, a11, a12, a21, a22, b11, b12, b21, b22, s1, s2, s3, s4, t1, t2, t3, t4, p1, p2, p3, p4, p5, p6, p7);
        // naive (A, B, D, n);
        // printf("C is:\n");
        // print_matrix (C, n);
    
        // printf("D is:\n");
        // print_matrix (D, n);

        // TODO: free memory

    }
    return 0;
}

void strassen (int size, double *A, double *B, double *C, double *a11, 
                  double *a12, double *a21, double *a22, double *b11, double *b12,
                  double *b21, double *b22, double *s1, double *s2, double *s3,
                  double *s4, double *t1, double *t2, double *t3, double *t4, 
                  double *p1, double *p2, double *p3, double *p4, double *p5, 
                  double *p6, double *p7)
{
    // tr == tc == 2. 
    if (size == 2) 
    {
        C[0] = (A[0] * B[0]) + (A[1] * B[2]);
        C[1] = (A[0] * B[1]) + (A[1] * B[3]);
        C[2] = (A[2] * B[0]) + (A[3] * B[2]);
        C[3] = (A[2] * B[1]) + (A[3] * B[3]);
        return;
    }
    int newSize = size / 2;

    for (int i = 0; i < newSize; i++) 
    {
        for (int j = 0; j < newSize; j++) 
        {
            a11[i * newSize + j] = A[i * size + j]; 
            a12[i * newSize + j] = A[i * size + (j + newSize)];
            a21[i * newSize + j] = A[(i + newSize) * size + j];
            a22[i * newSize + j] = A[(i + newSize) * size + (j + newSize)];

            b11[i * newSize + j] = B[i * size + j]; 
            b12[i * newSize + j] = B[i * size + (j + newSize)];
            b21[i * newSize + j] = B[(i + newSize) * size + j];
            b22[i * newSize + j] = B[(i + newSize) * size + (j + newSize)];
        }
    }
    // printf("A is: \n");
    // print_matrix(A, size);
    // printf("a11 is:\n");
    // print_matrix(a11, newSize);
    // printf("a12 is:\n");
    // print_matrix(a12, newSize);
    // printf("a21 is:\n");
    // print_matrix(a21, newSize);
    // printf("a22 is:\n");
    // print_matrix(a22, newSize);
    // printf("STOP\n");

    matrix_add   (a21, a22, s1, newSize);
    matrix_subtr (s1, a11,  s2, newSize);
    matrix_subtr (a11, a21, s3, newSize);
    matrix_subtr (a12, s2,  s4, newSize);

    matrix_subtr (b12, b11, t1, newSize);
    matrix_subtr (b22, t1,  t2, newSize);
    matrix_subtr (b22, b12, t3, newSize);
    matrix_subtr (b21, t2,  t4, newSize);

    strassen(newSize, a11, b11, p1, a11, a12, a21, a22, b11, b12, b21, b22, s1, s2, s3, s4, t1, t2, t3, t4, p1, p2, p3, p4, p5, p6, p7); 
    strassen(newSize, a12, b21, p2, a11, a12, a21, a22, b11, b12, b21, b22, s1, s2, s3, s4, t1, t2, t3, t4, p1, p2, p3, p4, p5, p6, p7);
    strassen(newSize, s1, t1, p3  , a11, a12, a21, a22, b11, b12, b21, b22, s1, s2, s3, s4, t1, t2, t3, t4, p1, p2, p3, p4, p5, p6, p7);
    strassen(newSize, s2, t2, p4  , a11, a12, a21, a22, b11, b12, b21, b22, s1, s2, s3, s4, t1, t2, t3, t4, p1, p2, p3, p4, p5, p6, p7); 
    strassen(newSize, s3, t3, p5  , a11, a12, a21, a22, b11, b12, b21, b22, s1, s2, s3, s4, t1, t2, t3, t4, p1, p2, p3, p4, p5, p6, p7);
    strassen(newSize, s4, b22, p6 , a11, a12, a21, a22, b11, b12, b21, b22, s1, s2, s3, s4, t1, t2, t3, t4, p1, p2, p3, p4, p5, p6, p7);
    strassen(newSize, a22, t4, p7 , a11, a12, a21, a22, b11, b12, b21, b22, s1, s2, s3, s4, t1, t2, t3, t4, p1, p2, p3, p4, p5, p6, p7);    

    // join partitions back into C
    for (int i = 0; i < newSize; i++) 
    {
        for (int j = 0; j < newSize; j++) 
        {
            C[i * size + j]                         = p1[i * newSize + j] + p2[i * newSize + j];
            C[i * size + (j + newSize)]             = p1[i * newSize + j] + p4[i * newSize + j] + p3[i * newSize + j] + p6[i * newSize + j]; 
            C[(i + newSize) * size + j]             = p1[i * newSize + j] + p4[i * newSize + j] + p5[i * newSize + j] + p7[i * newSize + j]; 
            C[(i + newSize) * size + (j + newSize)] = p1[i * newSize + j] + p4[i * newSize + j] + p5[i * newSize + j] + p3[i * newSize + j];
        }
    }
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
            printf("for (%d, %d), morton is %d\n", row, col, res);
            morton[res] = matrix[row * size + col];
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
    for (uint32_t i = 0; i < size; i++)
    {
        for (uint32_t j = 0; j < size; j++)
        {
            res[i * size + j] = a[i * size + j] - b[i * size + j];
        }
    }
}

void matrix_add (double *a, double *b, double *res, uint32_t size)
{
    for (uint32_t i = 0; i < size; i++)
    {
        for (uint32_t j = 0; j < size; j++)
        {
            res[i * size + j] = a[i * size + j] + b[i * size + j];
        }
    }
}

void naive (double *a, double *b, double *res, uint32_t size)
{
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
            res[i * size + j] += sum;
        }
    }
}

double rand_from (double min, double max)
{
    double rand_decimal = ((double) rand() / RAND_MAX) * (max - min) + min;
    
    return (uint32_t) (rand_decimal * 100.0 + 0.5) / 100.0;
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