#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>

void  strassen(int, double *, double *, double *, double *, 
              double *, double *,
              double *, double *, double *, double *, double *, 
              double *, double *, double *, double *, double *, 
              double *, double *);
void         morton_naive           (double *, double *, double *, uint32_t);
void         matrix_subtr    (double *, double *, double *, uint32_t);
void         matrix_add      (double *, double *, double *, uint32_t);
void         naive           (double *, double *, double *, uint32_t);
unsigned int S               (unsigned int, unsigned int);
void         convertToMorton (double *, double *, int);
void         print_matrix    (double *, uint32_t);
void         spcl_print      (double *, int, int);
double       rand_from       (double, double);
int          L_r             (int, int, int);  
void         print_bits      (unsigned int);
int          layout          (int, int);

const unsigned int E[] = {0x55555555, 0x33333333, 0x0F0F0F0F, 0x00FF00FF};
const unsigned int F[] = {1, 2, 4, 8};
const unsigned int TR = 2, TC = 2;
int sizes[] = {8};
int main()
{
    for (int s = 0; s < 1; s++)
    {
        int n = sizes[s];
        int newSize = n / 2;
        // index into A as A[i * cur_size + j]
        double *A   = (double *) calloc (n * n, sizeof (double));
        double *B   = (double *) calloc (n * n, sizeof (double));
        double *A_Morton =  (double *) calloc (n * n, sizeof (double));
        double *B_Morton =  (double *) calloc (n * n, sizeof (double));
        double *C   = (double *) calloc (n * n, sizeof (double));
        double *D   = (double *) calloc (n * n, sizeof (double));
        
        double *a11 = (double *) calloc (newSize * newSize, sizeof (double));
        double *a12 = (double *) calloc (newSize * newSize, sizeof (double)); 
        double *a21 = (double *) calloc (newSize * newSize, sizeof (double));
        double *a22 = (double *) calloc (newSize * newSize, sizeof (double));

        double *b11 = (double *) calloc (newSize * newSize, sizeof (double));
        double *b12 = (double *) calloc (newSize * newSize, sizeof (double)); 
        double *b21 = (double *) calloc (newSize * newSize, sizeof (double));
        double *b22 = (double *) calloc (newSize * newSize, sizeof (double));

        double *s1  = (double *) calloc (newSize * newSize, sizeof (double));
        double *s2  = (double *) calloc (newSize * newSize, sizeof (double));
        double *s3  = (double *) calloc (newSize * newSize, sizeof (double));
        double *s4  = (double *) calloc (newSize * newSize, sizeof (double));

        double *t1  = (double *) calloc (newSize * newSize, sizeof (double));
        double *t2  = (double *) calloc (newSize * newSize, sizeof (double));
        double *t3  = (double *) calloc (newSize * newSize, sizeof (double));
        double *t4  = (double *) calloc (newSize * newSize, sizeof (double));

        double *p1  = (double *) calloc (newSize * newSize, sizeof (double));
        double *p3  = (double *) calloc (newSize * newSize, sizeof (double));
        double *p2  = (double *) calloc (newSize * newSize, sizeof (double));
        double *p5  = (double *) calloc (newSize * newSize, sizeof (double));
        double *p4  = (double *) calloc (newSize * newSize, sizeof (double));
        double *p7  = (double *) calloc (newSize * newSize, sizeof (double));
        double *p6  = (double *) calloc (newSize * newSize, sizeof (double));


        // initialize A and B 
        for (int i = 0; i < n * n; i++) 
        {
            A[i] = i;
            B[i] = i;
            // for (int j = 0; j < n; j++) {
            //     A[i * n + j] = i + j;
            //     B[i * n + j] = i + j;
            // }
            
            // for (int j = 0; j < n; j++) 
            // {
            //     A[i * n + j] = i + j;
            //     B[i * n + j] = j;
            // }
        }

        convertToMorton(A, A_Morton, n);
        // printf("Initially A_Morton is:\n");
        // print_matrix (A_Morton, n);
        convertToMorton(B, B_Morton, n);
        // printf("Initially B_Morton is:\n");
        printf("A is \n");
        print_matrix(A, n);
        printf("B is \n");
        print_matrix(B, n);
        // print_matrix (B_Morton, n);
        strassen (n, A_Morton, B_Morton, C, s1, s2, s3, s4, t1, t2, t3, t4, p1, p2, p3, p4, p5, p6, p7);
        naive (A, B, D, n);
        printf("C is:\n");
        print_matrix (C, n);
    
        printf("D is:\n");
        print_matrix (D, n);

        // TODO: free memory

    }
    return 0;
}

void strassen (int size, double *A, double *B, double *C,
                  double *s1, double *s2, double *s3,
                  double *s4, double *t1, double *t2, double *t3, double *t4, 
                  double *p1, double *p2, double *p3, double *p4, double *p5, 
                  double *p6, double *p7)
{
    // tr == tc == 2. 
    if (size == TR) 
    {
        naive   (A, B, C, size);
        return;
    }
    int newSize = size / 2;
    printf("__A is__ \n");
    print_matrix(A, size);
    printf("__B__ is \n");
    print_matrix(B, size);
    
    // A[0:15]
    printf("**a11** is \n");
    spcl_print(A, 0, newSize);
    // A[16:31]
    printf("\n**a12** is \n");
    spcl_print(A, newSize * newSize, newSize);
    
    // A[32:47]
    printf("\n**a21** is \n");
    spcl_print(A, newSize * newSize * 2, newSize);
    
    //A[48:64]
    printf("\n**a22** is \n");
    spcl_print(A, newSize * newSize * 3, newSize);

    // a11 = A + 0
    // a12 = A + ((newSize * newSize) * 1)
    // a21 = A + ((newSize * newSize) * 2)
    // a22 = A + ((newSize * newSize) * 3)
    
    // S1 = A21 + A22
    matrix_add   (A + ((newSize * newSize) * 2) , A + ((newSize * newSize) * 3), s1, newSize);
    // S2 = S1 - A11
    matrix_subtr (s1, A,  s2, newSize);
    // S3 = A11 - A21
    matrix_subtr (A, A + ((newSize * newSize) * 2), s3, newSize);
    // S4 = A12 - S2
    matrix_subtr (A + ((newSize * newSize) * 1), s2,  s4, newSize);

    printf("\n**s1** is \n");
    print_matrix(s1, newSize);
    
    printf("**s2** is \n");
    print_matrix(s2, newSize);

    printf("**s3** is \n");
    print_matrix(s3, newSize);

    printf("**s4** is \n");
    print_matrix(s4, newSize);

    // T1 = B12 - B11
    matrix_subtr (B + ((newSize * newSize) * 1), B, t1, newSize);
    // T2 = B22 - T1
    matrix_subtr (B + ((newSize * newSize) * 3), t1,  t2, newSize);
    // T3 = B22 - B12
    matrix_subtr (B + ((newSize * newSize) * 3), B + ((newSize * newSize) * 1), t3, newSize);
    // T4 = B21 - T2
    matrix_subtr (B + ((newSize * newSize) * 2), t2,  t4, newSize);

    printf("**t1** is \n");
    print_matrix(t1, newSize);
    printf("**t2** is \n");
    print_matrix(t2, newSize);
    printf("**t3** is \n");
    print_matrix(t3, newSize);
    printf("**t4** is \n");
    print_matrix(t4, newSize);

    // a11 = A + 0
    // a12 = A + ((newSize * newSize) * 1)
    // a21 = A + ((newSize * newSize) * 2)
    // a22 = A + ((newSize * newSize) * 3)

    // P1 = A11 * B11
    strassen(newSize, A, B, p1, s1, s2, s3, s4, t1, t2, t3, t4, p1, p2, p3, p4, p5, p6, p7); 
    // P2 = A12 * B21
    strassen(newSize, A + ((newSize * newSize) * 1), B + ((newSize * newSize) * 2), p2, s1, s2, s3, s4, t1, t2, t3, t4, p1, p2, p3, p4, p5, p6, p7);
    // P3 = S1 * T1
    strassen(newSize, s1, t1, p3, s1, s2, s3, s4, t1, t2, t3, t4, p1, p2, p3, p4, p5, p6, p7);
    // P4 = S2 * T2
    strassen(newSize, s2, t2, p4, s1, s2, s3, s4, t1, t2, t3, t4, p1, p2, p3, p4, p5, p6, p7); 
    // P5 = S3 * T3
    strassen(newSize, s3, t3, p5, s1, s2, s3, s4, t1, t2, t3, t4, p1, p2, p3, p4, p5, p6, p7);
    // P6 = S4 * B22
    strassen(newSize, s4, B + ((newSize * newSize) * 3), p6, s1, s2, s3, s4, t1, t2, t3, t4, p1, p2, p3, p4, p5, p6, p7);
    // P7 = A22 * T4
    strassen(newSize, A + ((newSize * newSize) * 3), t4, p7, s1, s2, s3, s4, t1, t2, t3, t4, p1, p2, p3, p4, p5, p6, p7);    

    printf("**p1** is \n");
    print_matrix(p1, newSize);
    printf("**p2** is \n");
    print_matrix(p2, newSize);
    printf("**p3** is \n");
    print_matrix(p3, newSize);
    printf("**p4** is \n");
    print_matrix(p4, newSize);
    printf("**p5** is \n");
    print_matrix(p5, newSize);
    printf("**p6** is \n");
    print_matrix(p6, newSize);
    printf("**p7** is \n");
    print_matrix(p7, newSize);

    // a11 = A + 0
    // a12 = A + ((newSize * newSize) * 1)
    // a21 = A + ((newSize * newSize) * 2)
    // a22 = A + ((newSize * newSize) * 3)

    // join partitions back into C
    for (int i = 0; i < newSize * newSize; i++) {
        C[i + 0] = p1[i] + p2[i];
        C[i + ((newSize * newSize) * 1)] = p1[i] + p4[i] + p3[i] + p6[i];
        C[i + ((newSize * newSize) * 2)] = p1[i] + p4[i] + p5[i] + p7[i];
        C[i + ((newSize * newSize) * 3)] = p1[i] + p4[i] + p5[i] + p3[i];
    }

    // for (int i = 0; i < newSize; i++) 
    // {
    //     for (int j = 0; j < newSize; j++) 
    //     {
    //         C[i * newSize + j] = p1[i * newSize + j] + p2[i * newSize + j];
    //         C[i * newSize + (j + newSize)] = p1[i * newSize + j] + p4[i * newSize + j] + p3[i * newSize + j] + p6[i * newSize + j]; 
    //         C[(i + newSize) * newSize + j] = p1[i * newSize + j] + p4[i * newSize + j] + p5[i * newSize + j] + p7[i * newSize + j]; 
    //         C[(i + newSize) * newSize + (j + newSize) * newSize] = p1[i * newSize + j] + p4[i * newSize + j] + p5[i * newSize + j] + p3[i * newSize + j];
    //     }
    // }
    
    printf("at the end, C is \n");
    print_matrix(C, size);
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

void morton_naive (double *a, double *b, double *res, uint32_t size) {
    double sum = 0.0;
    for (uint32_t i = 0; i < size; i++)
    {
        for (uint32_t j = 0; j < size; j++)
        {
            sum = 0.0;
            for (uint32_t k = 0; k < size; k++)
            {
                sum += a[layout(i, k)] * b[layout(k, j)];
            }
            res[layout(i, j)] += sum;
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