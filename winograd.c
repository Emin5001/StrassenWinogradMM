#include <stdio.h>
#include <stdlib.h>
// #include <papi.h>
#include <time.h>


// int PAPI_one[3] = 
// {
//     PAPI_TOT_INS,
//     PAPI_LST_INS,
//     PAPI_FP_INS,
// };

// int PAPI_two[4] = 
// {
//     PAPI_L1_DCA,
//     PAPI_L1_DCM,
//     PAPI_L2_DCA,
//     PAPI_L2_DCM,
// };

// int PAPI_three[1] = { PAPI_TOT_CYC };

void strassen(int, double **, double **, double **, double **, double **, 
              double **, double **, double **, double **, double **, 
              double **, double **, double **, double **, double **,
              double **, double **, double **, double **, double **, 
              double **, double **, double **, double **, double **, 
              double **, double **, double **, double **, double **, 
              double **, double **, double **);
void naive(double **, double **, double **, int);
double rand_from (double, double);
void reset_matrices(double **, double **, double **, int);
void print_one_data (long long *, int);
void print_two_data (long long *, int);
void print_three_data (long long *, int);
void matrix_add (double **, double **, double **, int);
void matrix_subtract (double **, double **, double **, int);
void strassen_temp (double **, double **, double **, int);
void print_matrix (double **, int);

int sizes[12] = {2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096};
int main()
{
    int SIZE = 8;
    int half = SIZE / 2;

    // malloc all of the components we will need for strassen-winograd
    double **A  = malloc (SIZE * sizeof(double *));
    double **B  = malloc (SIZE * sizeof(double *));
    double **C  = malloc (SIZE * sizeof(double *));
    double **D  = malloc (SIZE * sizeof(double *)); // used for checking validity.

    double **p1 = malloc (half * sizeof(double *));
    double **p2 = malloc (half * sizeof(double *));
    double **p3 = malloc (half * sizeof(double *));
    double **p4 = malloc (half * sizeof(double *));
    double **p5 = malloc (half * sizeof(double *));
    double **p6 = malloc (half * sizeof(double *));
    double **p7 = malloc (half * sizeof(double *));

    double **u1 = malloc (half * sizeof(double *));
    double **u2 = malloc (half * sizeof(double *));
    double **u3 = malloc (half * sizeof(double *));
    double **u4 = malloc (half * sizeof(double *));
    double **u5 = malloc (half * sizeof(double *));
    double **u6 = malloc (half * sizeof(double *));
    double **u7 = malloc (half * sizeof(double *));

    double **s1 = malloc (half * sizeof(double *));
    double **s2 = malloc (half * sizeof(double *));
    double **s3 = malloc (half * sizeof(double *));
    double **s4 = malloc (half * sizeof(double *));

    double **t1 = malloc (half * sizeof(double *));
    double **t2 = malloc (half * sizeof(double *));
    double **t3 = malloc (half * sizeof(double *));
    double **t4 = malloc (half * sizeof(double *));

    double **a11 = malloc (half * sizeof(double *));
    double **a12 = malloc (half * sizeof(double *));
    double **a21 = malloc (half * sizeof(double *));
    double **a22 = malloc (half * sizeof(double *));

    double **b11 = malloc (half * sizeof(double *));
    double **b12 = malloc (half * sizeof(double *));
    double **b21 = malloc (half * sizeof(double *));
    double **b22 = malloc (half * sizeof(double *));

    for (int i = 0; i < SIZE; i++)
    {
        A[i]  = malloc(SIZE * sizeof(double));
        B[i]  = malloc(SIZE * sizeof(double));
        C[i]  = malloc(SIZE * sizeof(double));
        D[i]  = malloc(SIZE * sizeof(double));
        for (int j = 0; j < SIZE; j++)
        {
            A[i][j] = i;
            B[i][j] = j;
            C[i][j] = 0;
            D[i][j] = 0;
        }
    }

    for (int i = 0; i < half; i++)
    {
        p1[i] = malloc(half * sizeof(double));
        p2[i] = malloc(half * sizeof(double));
        p3[i] = malloc(half * sizeof(double));
        p4[i] = malloc(half * sizeof(double));
        p5[i] = malloc(half * sizeof(double));
        p6[i] = malloc(half * sizeof(double));
        p7[i] = malloc(half * sizeof(double));

        u1[i] = malloc(half * sizeof(double));
        u2[i] = malloc(half * sizeof(double));
        u3[i] = malloc(half * sizeof(double));
        u4[i] = malloc(half * sizeof(double));
        u5[i] = malloc(half * sizeof(double));
        u6[i] = malloc(half * sizeof(double));
        u7[i] = malloc(half * sizeof(double));

        s1[i] = malloc(half * sizeof(double));
        s2[i] = malloc(half * sizeof(double));
        s3[i] = malloc(half * sizeof(double));
        s4[i] = malloc(half * sizeof(double));

        t1[i] = malloc(half * sizeof(double));
        t2[i] = malloc(half * sizeof(double));
        t3[i] = malloc(half * sizeof(double));
        t4[i] = malloc(half * sizeof(double));

        a11[i] = malloc(half * sizeof(double));
        a12[i] = malloc(half * sizeof(double));
        a21[i] = malloc(half * sizeof(double));
        a22[i] = malloc(half * sizeof(double));

        b11[i] = malloc(half * sizeof(double));
        b12[i] = malloc(half * sizeof(double));
        b21[i] = malloc(half * sizeof(double));
        b22[i] = malloc(half * sizeof(double));
    }

    // init all components that we will need to 0
    for (int i = 0; i < half; i++)
    {
        for (int j = 0; j < half; j++)
        {
            p1[i][j]  = 0;
            p2[i][j]  = 0;
            p3[i][j]  = 0;
            p4[i][j]  = 0;
            p5[i][j]  = 0;
            p6[i][j]  = 0;
            p7[i][j]  = 0;

            a11[i][j] = 0;
            a12[i][j] = 0;
            a21[i][j] = 0;
            a22[i][j] = 0;

            b11[i][j] = 0;
            b12[i][j] = 0;
            b21[i][j] = 0;
            b22[i][j] = 0;

             s1[i][j] = 0;
             s2[i][j] = 0;
             s3[i][j] = 0;
             s4[i][j] = 0;

             t1[i][j] = 0;
             t2[i][j] = 0;
             t3[i][j] = 0;
             t4[i][j] = 0;

             u1[i][j] = 0;
             u2[i][j] = 0;
             u3[i][j] = 0;
             u4[i][j] = 0;
             u5[i][j] = 0;
             u6[i][j] = 0;
             u7[i][j] = 0;
        } 
    }
    printf("A is: \n");
    print_matrix (A, SIZE);
    
    printf("B is: \n");
    print_matrix (B, SIZE);
    // call strassen
    strassen(SIZE, A, B, C, a11, a12, a21, a22, b11, b12, b21, b22, s1, s2, s3, s4, t1, t2, t3, t4, p1, p2, p3, p4, p5, p6, p7, u1, u2, u3, u4, u5, u6, u7);
    // call naive; will use to compare.
    naive(A, B, D, SIZE);
    
    printf("\nC is: \n");
    print_matrix (C, SIZE);

    printf("D is: \n");
    print_matrix (D, SIZE);


    // for (int i = 0; i < SIZE; i++)
    // {
    //     for (int j = 0; j < SIZE; j++)
    //     {
    //         if (abs(C[i][j] - D[i][j]) != 0) 
    //         {
    //             printf("\033[1;31mError. C[%d][%d] is %f, \
    //             while D[%d][%d] is %.10f. \033[0m\n", i, j, C[i][j], i, j, D[i][j]);
    //         }
    //     }
    // }
    return 0;
}

// int main ()
// {
    
//     int size = 4;
//     double **A = malloc (size * sizeof (double *));
//     double **B = malloc (size * sizeof (double *));
//     double **C = malloc (size * sizeof (double *));
//     double **D = malloc (size * sizeof (double *));

//     for (int i = 0; i < size; i++)
//     {
//         A[i] = malloc (size * sizeof (double));
//         B[i] = malloc (size * sizeof (double));
//         C[i] = malloc (size * sizeof (double));
//         D[i] = malloc (size * sizeof (double));
        
//         for (int j = 0; j < size; j++)
//         {
//             A[i][j] = rand_from (1, 3);
//             B[i][j] = rand_from (1, 3);
//             C[i][j] = 0;
//             D[i][j] = 0;
//         }
//     }

//     strassen_temp (A, B, C, size);
//     naive (A, B, D, size);
//     for (int i = 0; i < size; i++)
//     {
//         for (int j = 0; j < size; j++)
//         {
//             if (abs(C[i][j] - D[i][j]) != 0) 
//             {
//                 printf("\033[1;31mError. C[%d][%d] is %f, \
//                 while D[%d][%d] is %.10f. \033[0m\n", i, j, C[i][j], i, j, D[i][j]);
//             }
//         }
//     }


//     return 0;
// }

void strassen_temp(double **A, double **B, double **C, int size)
{
    if (size == 2)
    {
        C[0][0] = (A[0][0] * B[0][0]) + (A[0][1] * B[1][0]);
        C[0][1] = (A[0][0] * B[0][1]) + (A[0][1] * B[1][1]);
        C[1][0] = (A[1][0] * B[0][0]) + (A[1][1] * B[1][0]);
        C[1][1] = (A[1][0] * B[0][1]) + (A[1][1] * B[1][1]);
    } else
    {
        int half = size / 2;

        double **a11 = malloc (half * sizeof(double *));
        double **a12 = malloc (half * sizeof(double *));
        double **a21 = malloc (half * sizeof(double *));
        double **a22 = malloc (half * sizeof(double *));

        double **b11 = malloc (half * sizeof(double *));
        double **b12 = malloc (half * sizeof(double *));
        double **b21 = malloc (half * sizeof(double *));
        double **b22 = malloc (half * sizeof(double *));

        double **c11 = malloc (half * sizeof(double *));
        double **c12 = malloc (half * sizeof(double *));
        double **c21 = malloc (half * sizeof(double *));
        double **c22 = malloc (half * sizeof(double *));

        double **s1 = malloc (half * sizeof(double *));
        double **s2 = malloc (half * sizeof(double *));
        double **s3 = malloc (half * sizeof(double *));
        double **s4 = malloc (half * sizeof(double *));

        double **t1 = malloc (half * sizeof(double *));
        double **t2 = malloc (half * sizeof(double *));
        double **t3 = malloc (half * sizeof(double *));
        double **t4 = malloc (half * sizeof(double *));

        double **p1 = malloc (half * sizeof(double *));
        double **p2 = malloc (half * sizeof(double *));
        double **p3 = malloc (half * sizeof(double *));
        double **p4 = malloc (half * sizeof(double *));
        double **p5 = malloc (half * sizeof(double *));
        double **p6 = malloc (half * sizeof(double *));
        double **p7 = malloc (half * sizeof(double *));

        double **u1 = malloc (half * sizeof(double *));
        double **u2 = malloc (half * sizeof(double *));
        double **u3 = malloc (half * sizeof(double *));
        double **u4 = malloc (half * sizeof(double *));
        double **u5 = malloc (half * sizeof(double *));
        double **u6 = malloc (half * sizeof(double *));
        double **u7 = malloc (half * sizeof(double *));

        for (int i = 0; i < half; i++)
        {
            p1[i] = malloc(half * sizeof(double));
            p2[i] = malloc(half * sizeof(double));
            p3[i] = malloc(half * sizeof(double));
            p4[i] = malloc(half * sizeof(double));
            p5[i] = malloc(half * sizeof(double));
            p6[i] = malloc(half * sizeof(double));
            p7[i] = malloc(half * sizeof(double));

            u1[i] = malloc(half * sizeof(double));
            u2[i] = malloc(half * sizeof(double));
            u3[i] = malloc(half * sizeof(double));
            u4[i] = malloc(half * sizeof(double));
            u5[i] = malloc(half * sizeof(double));
            u6[i] = malloc(half * sizeof(double));
            u7[i] = malloc(half * sizeof(double));

            s1[i] = malloc(half * sizeof(double));
            s2[i] = malloc(half * sizeof(double));
            s3[i] = malloc(half * sizeof(double));
            s4[i] = malloc(half * sizeof(double));

            t1[i] = malloc(half * sizeof(double));
            t2[i] = malloc(half * sizeof(double));
            t3[i] = malloc(half * sizeof(double));
            t4[i] = malloc(half * sizeof(double));

            a11[i] = malloc(half * sizeof(double));
            a12[i] = malloc(half * sizeof(double));
            a21[i] = malloc(half * sizeof(double));
            a22[i] = malloc(half * sizeof(double));

            b11[i] = malloc(half * sizeof(double));
            b12[i] = malloc(half * sizeof(double));
            b21[i] = malloc(half * sizeof(double));
            b22[i] = malloc(half * sizeof(double));

            c11[i] = malloc(half * sizeof(double));
            c12[i] = malloc(half * sizeof(double));
            c21[i] = malloc(half * sizeof(double));
            c22[i] = malloc(half * sizeof(double));

            for (int j = 0; j < half; j++)
            {  
                a11[i][j] = 0;
                a12[i][j] = 0;
                a21[i][j] = 0;
                a22[i][j] = 0;

                b11[i][j] = 0;
                b12[i][j] = 0;
                b21[i][j] = 0;
                b22[i][j] = 0;

                c11[i][j] = 0;
                c12[i][j] = 0;
                c21[i][j] = 0;
                c22[i][j] = 0;
            }
        }

        for (int i = 0; i < half; i++)
        {
            for (int j = 0; j < half; j++)
            {
                a11[i][j] = A[i][j];
                a12[i][j] = A[i][j + half];
                a21[i][j] = A[i + half][j];
                a22[i][j] = A[i + half][j + half];

                b11[i][j] = B[i][j];
                b12[i][j] = B[i][j + half];
                b21[i][j] = B[i + half][j];
                b22[i][j] = B[i + half][j + half];

            }
        }

        matrix_add (s1, a21, a22, half);
        matrix_subtract (s2, s1, a11, half);
        matrix_subtract (s3, a11, a21, half);
        matrix_subtract (s4, a12, s2, half);

        matrix_subtract (t1, b12, b11, half);
        matrix_subtract (t2, b22, t1, half);
        matrix_subtract (t3, b22, b12, half);
        matrix_subtract (t4, b21, t2, half);

        strassen_temp (a11, b11, p1, half);
        strassen_temp (a12, b21, p2, half);
        strassen_temp (s1,  t1,   p3, half);
        strassen_temp (s2,  t2,   p4, half);
        strassen_temp (s3,  t3,  p5, half);
        strassen_temp (s4,  b22, p6, half);
        strassen_temp (a22, t4,  p7, half);

        matrix_add (u1, p1, p2, half);
        matrix_add (u2, p1, p4, half);
        matrix_add (u3, u2, p5, half);
        matrix_add (u4, u3, p7, half);
        matrix_add (u5, u3, p3, half);
        matrix_add (u6, u2, p3, half);
        matrix_add (u7, u6, p6, half);

        matrix_add (c11, p1, p2, half);
        matrix_add (c12, u3, p7, half);
        matrix_add (c22, u3, p3, half);
        matrix_add (c21, u6, p6, half);

        for (int i = 0; i < half; i++)
        {
            for (int j = 0; j < half; j++)
            {
                C[i][j] = c11[i][j];
                C[i][j + half] = c12[i][j];
                C[i + half][j] = c21[i][j];
                C[i + half][j + half] = c22[i][j];
            }
        }
    }   
}

void matrix_add (double **dest, double **a, double **b, int size)
{
    for (int r = 0; r < size; r++)
    {
        for (int c = 0; c < size; c++)
        {
            dest[r][c] = a[r][c] + b[r][c];
        }
    }
}

void matrix_subtract (double **dest, double **a, double **b, int size)
{
    for (int r = 0; r < size; r++)
    {
        for (int c = 0; c < size; c++)
        {
            dest[r][c] = a[r][c] - b[r][c];
        }
    }
}

void strassen(int size, double **A, double **B, double **C, double **a11, 
              double **a12, double **a21, double **a22, double **b11, 
              double **b12, double **b21, double **b22, double **s1, 
              double **s2, double **s3, double **s4, double **t1, 
              double **t2, double **t3, double **t4, double **p1, 
              double **p2, double **p3, double **p4, double **p5, 
              double **p6, double **p7, double **u1, double **u2, double **u3,
              double **u4, double **u5, double **u6, double **u7)
{
    if (size == 2)
    {
        printf("The address of the resulting array in the if statement is %p\n", (void *) C);
        C[0][0] = (A[0][0] * B[0][0]) + (A[0][1] * B[1][0]);
        C[0][1] = (A[0][0] * B[0][1]) + (A[0][1] * B[1][1]);
        C[1][0] = (A[1][0] * B[0][0]) + (A[1][1] * B[1][0]);
        C[1][1] = (A[1][0] * B[0][1]) + (A[1][1] * B[1][1]);
    } else
    {
        int half = size / 2;
        
        // calculate a11 --> a22 and b11 --> b22.
        for (int i = 0; i < half; i++)
        {
            for (int j = 0; j < half; j++)
            {
                a11[i][j] = A[i][j];
                b11[i][j] = B[i][j];

                a12[i][j] = A[i][j + half];
                b12[i][j] = B[i][j + half];

                a21[i][j] = A[i + half][j];
                b21[i][j] = B[i + half][j];

                a22[i][j] = A[i + half][j + half];
                b22[i][j] = B[i + half][j + half];
            }
        }
        printf("a11 is: \n");
        print_matrix(a11, half);
        
        printf("a12 is: \n");
        print_matrix(a12, half);
        
        printf("a21 is: \n");
        print_matrix(a21, half);
        
        printf("a22 is: \n");
        print_matrix(a22, half);

        printf("b11 is: \n");
        print_matrix(b11, half);
        
        printf("b12 is: \n");
        print_matrix(b12, half);
        
        printf("b21 is: \n");
        print_matrix(b21, half);
        
        printf("b22 is: \n");
        print_matrix(b22, half);
        
        // calculate s1 --> s4 and t1 --> t4
        matrix_add      (s1, a21, a22, half);
        matrix_subtract (s2, s1, a11, half);
        matrix_subtract (s3, a11, a21, half);
        matrix_subtract (s4, a12, s2, half);
        printf("\nS1 is: \n");
        print_matrix(s1, half);
        printf("S2 is: \n");
        print_matrix(s2, half);
        printf("S3 is: \n");
        print_matrix(s3, half);
        printf("S4 is: \n");
        print_matrix(s4, half);

        matrix_subtract (t1, b12, b11, half);
        matrix_subtract (t2, b22, t1, half);
        matrix_subtract (t3, b22, b12, half);
        matrix_subtract (t4, b21, t2, half);
        
        printf("\nT1 is: \n");
        print_matrix(t1, half);
        printf("T2 is: \n");
        print_matrix(t2, half);
        printf("T3 is: \n");
        print_matrix(t3, half);
        printf("T4 is: \n");
        print_matrix(t4, half);

        strassen(half, a11, b11,  p1, a11, a12, a21, a22, b11, b12, b21, b22, s1, s2, s3, s4, t1, t2, t3, t4, p1, p2, p3, p4, p5, p6, p7, u1, u2, u3, u4, u5, u6, u7); // p1 = a11 * b11
        strassen(half, a12, b21,  p2, a11, a12, a21, a22, b11, b12, b21, b22, s1, s2, s3, s4, t1, t2, t3, t4, p1, p2, p3, p4, p5, p6, p7, u1, u2, u3, u4, u5, u6, u7); // p2 = a12 * b21
        strassen(half, s1,  t1,  p3, a11, a12, a21, a22, b11, b12, b21, b22, s1, s2, s3, s4, t1, t2, t3, t4, p1, p2, p3, p4, p5, p6, p7, u1, u2, u3, u4, u5, u6, u7); // p3 = s1  * t1
        strassen(half, s2, t2,  p4, a11, a12, a21, a22, b11, b12, b21, b22, s1, s2, s3, s4, t1, t2, t3, t4, p1, p2, p3, p4, p5, p6, p7, u1, u2, u3, u4, u5, u6, u7); // p4 = s2  * t2
        strassen(half, s3,  t3,  p5, a11, a12, a21, a22, b11, b12, b21, b22, s1, s2, s3, s4, t1, t2, t3, t4, p1, p2, p3, p4, p5, p6, p7, u1, u2, u3, u4, u5, u6, u7); // p5 = s3  * t3
        strassen(half, s4,  b22, p6, a11, a12, a21, a22, b11, b12, b21, b22, s1, s2, s3, s4, t1, t2, t3, t4, p1, p2, p3, p4, p5, p6, p7, u1, u2, u3, u4, u5, u6, u7); // p6 = s4  * b22
        strassen(half, a22,  t4,  p7, a11, a12, a21, a22, b11, b12, b21, b22, s1, s2, s3, s4, t1, t2, t3, t4, p1, p2, p3, p4, p5, p6, p7, u1, u2, u3, u4, u5, u6, u7); // p7 = a22 * t4
        
        printf("\np1 is: \n");
        print_matrix (p1, half);
        printf("\np2 is: \n");
        print_matrix (p2, half);
        printf("\np3 is: \n");
        print_matrix (p3, half);
        printf("\np4 is: \n");
        print_matrix (p4, half);
        printf("\np5 is: \n");
        print_matrix (p5, half);
        printf("\np6 is: \n");
        print_matrix (p6, half);
        printf("\np7 is: \n");
        print_matrix (p7, half);

        // calculate u1 --> u7 & set C11 --> C22
        for (int i = 0; i < half; i++) 
        {
            for (int j = 0; j < half; j++) 
            {
                C[i][j]  = p1[i][j] + p2[i][j];                 // C11 = U1 = P1 + P2
                u2[i][j] = p1[i][j] + p4[i][j];
                u3[i][j] = u2[i][j] + p5[i][j];
                C[i + half][j] = u3[i][j] + p7[i][j];           // C21 = U4 = U3 + P7
                C[i + half][j + half] = u3[i][j] + p3[i][j];    // C22 = U5 = U3 + P3
                u6[i][j] = u2[i][j] + p3[i][j];
                C[i][j + half] = u6[i][j] + p6[i][j];           // C12 = U7 = U6 + P6
            }
        }
        printf("Now, C is:\n");
        print_matrix (C, size);
    }
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

double rand_from (double min, double max)
{
    double rand_decimal = ((double) rand() / RAND_MAX) * (max - min) + min;
    
    // rounds to nearest 2 decimal points.
    return (int)(rand_decimal * 100.0 + 0.5) / 100.0;
}

void reset_matrices(double **a, double **b, double **c, int size)
{
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            a[i][j] = rand_from(1, 3);
            b[i][j] = rand_from(1, 3);
            c[i][j] = 0.0;
        }
    }
}

void print_matrix (double **matrix, int size)
{
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            printf("%.4f ", matrix[i][j]);
        }
        printf("\n");
    }
}

void print_one_data (long long *counters, int size) 
{
    char names[3][20] = {"PAPI_TOT_INS", "PAPI_LST_INS", "PAPI_FP_INS"};
    printf("\n%-20s %-10s\n", "Name", "Value");
    printf("--------------------------------------------------\n");

    for (int i = 0; i < 3; i++) {
        printf("%-20s %-10lld\n", names[i], counters[i]);
    }

    printf("Now it is normalized data.\n");
    printf("\n%-20s %-10s\n", "Name", "Value");
    printf("--------------------------------------------------\n");

    for (int i = 0; i < 3; i++) 
    {
        printf("%-20s %.4f\n", names[i], ((double) counters[i] / (double) (size * size * size)));
    }
}

void print_two_data (long long *counters, int size) 
{
    printf("L1 Miss Rate is: %0.3f\n", ((float) counters[1] / (float) counters[0]));
    printf("L2 Miss Rate is: %0.3f\n", ((float) counters[3] / (float) counters[2]));
}

void print_three_data(long long *counters, int size) 
{
    char *name = "PAPI_TOT_CYC";
    printf("%-20s %-10lld\n\n", name, counters[0]);
}
