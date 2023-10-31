#include <stdio.h>
#include <time.h>
#include <stdlib.h>

// double** strassen (int, double **, double **, double **, double **, double **, double **, double **, double **, double **, double **, double **, double **, double **, double **, double **);

double* strassen (int size, double *A, double *B, double *C, double *a11, 
                  double *a12, double *a21, double *a22, double *b11, double *b12,
                  double *b21, double *b22, double *s1, double *s2, double *s3,
                  double *s4, double *t1, double *t2, double *t3, double *t4, 
                  double *p1, double *p2, double *p3, double *p4, double *p5, 
                  double *p6, double *p7);
    
double*  matrix_subtr (double *A, double *B, double *res, int size);
double*  matrix_add (double *A, double *B, double *res, int size);
void     naive (double *, double *, double *, int);
double   rand_from (double, double);
void     print_matrix (double *, int);

int main() {
    int sizes[] = {8};
    for (int s = 0; s < 1; s++)
    {
        int n = sizes[s];
        int newSize = n / 2;
        // index into A as A[i * cur_size + j]
        double *A   = (double *) malloc (n * n * sizeof (double));
        double *B   = (double *) malloc (n * n * sizeof (double));
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
        for (int i = 0; i < n; i++) 
        {
            for (int j = 0; j < n; j++) 
            {
                A[i * n + j] = i;
                B[i * n + j] = j;
            }
        }

        printf("Initially A is:\n");
        print_matrix (A, n);

        printf("Initially B is:\n");
        print_matrix (B, n);
        strassen (n, A, B, C, a11, a12, a21, a22, b11, b12, b21, b22, s1, s2, s3, s4, t1, t2, t3, t4, p1, p2, p3, p4, p5, p6, p7);
        naive (A, B, D, n);
        printf("C is:\n");
        print_matrix (C, n);
    
        printf("D is:\n");
        print_matrix (D, n);

        // TODO: free memory

    }

    return 0;
}
    


double* strassen(int size, double *A, double *B, double *C, double *a11, 
                 double *a12, double *a21, double *a22, double *b11, double *b12,
                 double *b21, double *b22, double *s1, double *s2, double *s3, 
                 double *s4, double *t1, double *t2, double *t3, double *t4, 
                 double *p1, double *p2, double *p3, double *p4, double *p5,
                 double *p6, double *p7) 
{
    if (size == 2) {
        C[0] = (A[0] * B[0]) + (A[1] * B[2]);
        C[1] = (A[0] * B[1]) + (A[1] * B[3]);
        C[2] = (A[2] * B[0]) + (A[3] * B[2]);
        C[3] = (A[2] * B[1]) + (A[3] * B[3]);
        return C;
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
    printf("A is: \n");
    print_matrix(A, size);
    printf("a11 is:\n");
    print_matrix(a11, newSize);
    printf("a12 is:\n");
    print_matrix(a12, newSize);
    printf("a21 is:\n");
    print_matrix(a21, newSize);
    printf("a22 is:\n");
    print_matrix(a22, newSize);
    printf("STOP\n");
    // for (int i = 0; i < newSize; i++)
    // {
    //     for (int j = 0; j < newSize; j++)
    //     {
    //         a11[i][j] = A[i][j];
    //     }
    // }

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

    return C;
}


void naive (double *A, double *B, double *res, int size)
{
    double sum = 0.0;
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            sum = 0.0;
            for (int k = 0; k < size; k++)
            {
                sum += A[i * size + k] * B[k * size + j];
            }
            res[i * size + j] += sum;
        }
    }
}

double* matrix_add (double *A, double *B, double *res, int size)
{
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            res[i * size + j] = A[i * size + j] + B[i * size + j];
        }
    }

    return res;
}

double* matrix_subtr (double *A, double *B, double *res, int size)
{
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            res[i * size + j] = A[i * size + j] - B[i * size + j];
        }
    }
    return res;
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

double rand_from (double min, double max)
{
    double rand_decimal = ((double) rand() / RAND_MAX) * (max - min) + min;
    
    // rounds to nearest 2 decimal points.
    return (int) (rand_decimal * 100.0 + 0.5) / 100.0;
}