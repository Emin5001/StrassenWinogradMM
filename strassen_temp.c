#include <stdio.h>
#include <stdlib.h>

// double** strassen (int, double **, double **, double **, double **, double **, double **, double **, double **, double **, double **, double **, double **, double **, double **, double **);

double** strassen (int, double **, double **, double **, double **, double **, 
                    double **, double **);

void naive (double **, double **, double **, int);
double** matrix_add (double **, double **, double **, int);
double** matrix_subtr (double **, double **, double **, int);
void print_matrix (double **, int);

int main() {
    int n = 8;
    int newSize = n / 2;

    double **A = (double **) malloc (n * sizeof (double *));
    double **B = (double **) malloc (n * sizeof (double *));
    double **C = (double **) malloc (n * sizeof (double *));
    double **D = (double **) malloc (n * sizeof (double *));
    
    // double **a11 = (double **) malloc (newSize * sizeof (double *));
    // double **a12 = (double **) malloc (newSize * sizeof (double *)); 
    // double **a21 = (double **) malloc (newSize * sizeof (double *));
    // double **a22 = (double **) malloc (newSize * sizeof (double *));

    // double **b11 = (double**)malloc(newSize*sizeof(double*));
    // double **b12 = (double**)malloc(newSize*sizeof(double*)); 
    // double **b21 = (double**)malloc(newSize*sizeof(double*));
    // double **b22 = (double**)malloc(newSize*sizeof(double*));

    double **c11 = (double **) malloc (newSize * sizeof (double *));
    double **c12 = (double **) malloc (newSize * sizeof (double *)); 
    double **c21 = (double **) malloc (newSize * sizeof (double *));
    double **c22 = (double **) malloc (newSize * sizeof (double *));

    // double **s1 = (double **) malloc (newSize * sizeof (double *));
    // double **s2 = (double **) malloc (newSize * sizeof (double *));
    // double **s3 = (double **) malloc (newSize * sizeof (double *));
    // double **s4 = (double **) malloc (newSize * sizeof (double *));

    // double **t1 = (double **) malloc (newSize * sizeof (double *));
    // double **t2 = (double **) malloc (newSize * sizeof (double *));
    // double **t3 = (double **) malloc (newSize * sizeof (double *));
    // double **t4 = (double **) malloc (newSize * sizeof (double *));

    // double **p1 = (double **) malloc (newSize * sizeof (double *));
    // double **p2 = (double **) malloc (newSize * sizeof (double *));
    // double **p3 = (double **) malloc (newSize * sizeof (double *));
    // double **p4 = (double **) malloc (newSize * sizeof (double *));
    // double **p5 = (double **) malloc (newSize * sizeof (double *));
    // double **p6 = (double **) malloc (newSize * sizeof (double *));
    // double **p7 = (double **) malloc (newSize * sizeof (double *));
    
    
    for (int i = 0; i < n; i++) 
    {
        A[i] = (double*) malloc (n * sizeof(double));
        B[i] = (double*) malloc (n * sizeof(double));
        C[i] = (double*) malloc (n * sizeof(double));
        D[i] = (double*) malloc (n * sizeof(double));
    }

    for (int i = 0; i < n / 2; i++)
    {
        // a11[i] = (double *) malloc (newSize * sizeof (double));
        // a12[i] = (double *) malloc (newSize * sizeof (double));
        // a21[i] = (double *) malloc (newSize * sizeof (double));
        // a22[i] = (double *) malloc (newSize * sizeof (double));

        // b11[i] = (double *) malloc (newSize * sizeof (double));
        // b12[i] = (double *) malloc (newSize * sizeof (double));
        // b21[i] = (double *) malloc (newSize * sizeof (double));
        // b22[i] = (double *) malloc (newSize * sizeof (double));

        c11[i] = (double *) malloc (newSize * sizeof (double));
        c12[i] = (double *) malloc (newSize * sizeof (double));
        c21[i] = (double *) malloc (newSize * sizeof (double));
        c22[i] = (double *) malloc (newSize * sizeof (double));

        // s1[i] = (double *) malloc (newSize * sizeof (double));
        // s2[i] = (double *) malloc (newSize * sizeof (double));
        // s3[i] = (double *) malloc (newSize * sizeof (double));
        // s4[i] = (double *) malloc (newSize * sizeof (double));

        // t1[i] = (double *) malloc (newSize * sizeof (double));
        // t2[i] = (double *) malloc (newSize * sizeof (double));
        // t3[i] = (double *) malloc (newSize * sizeof (double));
        // t4[i] = (double *) malloc (newSize * sizeof (double));

        // p1[i] = (double *) malloc (newSize * sizeof (double));
        // p2[i] = (double *) malloc (newSize * sizeof (double));
        // p3[i] = (double *) malloc (newSize * sizeof (double));
        // p4[i] = (double *) malloc (newSize * sizeof (double));
        // p5[i] = (double *) malloc (newSize * sizeof (double));
        // p6[i] = (double *) malloc (newSize * sizeof (double));
        // p7[i] = (double *) malloc (newSize * sizeof (double));
    }

    // initialize A and B 
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            A[i][j] = i;
            B[i][j] = j;
        }
    }
    printf("Initially A is:\n");
    print_matrix (A, n);

    printf("Initially B is:\n");
    print_matrix (B, n);
    // strassen (n, A, B, C, s1, s2, s3, s4, t1, t2, t3, t4, c11, c12, c21, c22);
    strassen (n, A, B, C, c11, c12, c21, c22);

    naive (A, B, D, n);
    // print C
    printf("C is:\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%lf ", C[i][j]); 
        }
        printf("\n");
    }
  
    printf("D is:\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%lf ", D[i][j]); 
        }
        printf("\n");
    }

    // free memory
    for (int i = 0; i < n; i++) {
        free(A[i]);
        free(B[i]);
        free(C[i]);
    }
    free(A);
    free(B);
    free(C);

    return 0;
}


double** strassen(int n, double **A, double **B, double **C, double **c11, 
                  double **c12, double **c21, double **c22) {
    if (n == 2) {
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

    // double **c11 = (double **) malloc (newSize * sizeof (double *));
    // double **c12 = (double **) malloc (newSize * sizeof (double *));
    // double **c21 = (double **) malloc (newSize * sizeof (double *));
    // double **c22 = (double **) malloc (newSize * sizeof (double *));

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

        // c11[i] = (double *) malloc (newSize * sizeof (double));
        // c12[i] = (double *) malloc (newSize * sizeof (double));
        // c21[i] = (double *) malloc (newSize * sizeof (double));
        // c22[i] = (double *) malloc (newSize * sizeof (double));

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
    printf("\n\n");
    // divide matrices
    printf("\nnewSize is %d", newSize);
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
    
    printf("\nA is:\n");
    print_matrix (A, n);
    for (int i = 0; i < newSize; i++)
    {
        for (int j = 0; j < newSize; j++)
        {
            a11[i][j] = A[i][j];
        }
    }
    printf("a11 is:\n");
    print_matrix (a11, newSize);

    printf("a12 is:\n");
    print_matrix (a12, newSize);
    
    printf("a21 is:\n");
    print_matrix (a21, newSize);

    printf("a22 is:\n");
    print_matrix (a22, newSize);

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


void naive (double **A, double **B, double **res, int size)
{
    double sum = 0.0;
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            sum = 0.0;
            for (int k = 0; k < size; k++)
            {
                sum += A[i][k] * B[k][j];
            }
            res[i][j] += sum;
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