#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void  strassen(int, double *, double *, double *, double *, 
              double *, double *, double *, double *, double *,
              double *, double *, double *, double *, double *,
              double *, double *, double *, double *, double *, 
              double *, double *, double *, double *, double *, 
              double *p6, double *p7);

void   matrix_subtr    (double *, double *, double *, int);
void   matrix_add      (double *, double *, double *, int);
void   naive           (double *, double *, double *, int);
void   convertToMorton (double *, double *, int);
double rand_from       (double, double);
void   print_matrix    (double *, int);
long long    interleave_bits (int, int, int);
void   print_bits      (unsigned int);
int main()
{
    int size = 8;
    double *A = (double *) malloc (size * size * sizeof (double));
    double *morton = (double *) malloc (size * size * sizeof (double));
    for (int i = 0; i < size * size; i++)
    {
        A[i] = i;
    }    

    printf("Initially the matrix is: \n");
    print_matrix(A, size);
    convertToMorton(A, morton, size);
    printf("Now morton is: \n");
    print_matrix(morton, size);

    return 0;
}


void convertToMorton(double *matrix, double *morton, uint_fast32_t size)
{
    for (uint_fast32_t row = 0; row < size; row++)
    {
        for (uint_fast32_t col = 0; col < size; col++)
        {
            // int mortonIndex = interleave_bits(row, col, size);
            // printf("row: %d col: %d morton: %d\n", row, col, mortonIndex);
            // morton[mortonIndex] = matrix[row * size + col];
        }
    }
}

// r=0, c=2 
// 00000000
// 00000010 ==> 0000000000000100 => 4. however, if size = 8, it should be
long long interleave_bits (uint_fast32_t x, uint_fast32_t y, int size)
{
    long long z = 0;
    for (int i = 0; i < sizeof(x) * 8 / log(size); i++)
    {
        z |= (x & 1U << i) << i | (y & 1U << i) << (i + 1);
    }

    return z;
}

void matrix_subtr (double *a, double *b, double *res, int size)
{
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            res[i * size + j] = a[i * size + j] - b[i * size + j];
        }
    }
}

void matrix_add (double *a, double *b, double *res, int size)
{
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            res[i * size + j] = a[i * size + j] + b[i * size + j];
        }
    }
}

void naive (double *a, double *b, double *res, int size)
{
    double sum = 0.0;
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            sum = 0.0;
            for (int k = 0; k < size; k++)
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
    
    return (int) (rand_decimal * 100.0 + 0.5) / 100.0;
}

void print_matrix (double *matrix, int size)
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

// 3 --> 00000000 00000000 00000000 00000011 
//     & 10000000 00000000 00000000 00000000 
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