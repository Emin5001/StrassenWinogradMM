#include <stdio.h>
#include <stdlib.h>
#include <papi.h>
#include <time.h>


int PAPI_one[3] = 
{
    PAPI_TOT_INS,
    PAPI_LST_INS,
    PAPI_FP_INS,
};

int PAPI_two[4] = 
{
    PAPI_L1_DCA,
    PAPI_L1_DCM,
    PAPI_L2_DCA,
    PAPI_L2_DCM,
};

int PAPI_three[1] = { PAPI_TOT_CYC };

void strassen(int, double **, double **, double **, double **, double **, 
              double **, double **, double **, double **, double **, 
              double **, double **, double **, double **, double **,
              double **, 
              double **, double **, double **, double **, 
              double **, double **, double **, double **, 
              double **, double **, double **, double **, 
              double **, double **, double **, double **, double **);
void naive(double **, double **, double **, int);
double rand_from (double, double);
void reset_matrices(double **, double **, double **, int);
void print_one_data (long long *, int);
void print_two_data (long long *, int);
void print_three_data (long long *, int);

int sizes[12] = {2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096};
int main()
{
    if (PAPI_library_init(PAPI_VER_CURRENT) < 0) 
    {
        fprintf(stderr, "Error intializing PAPI!");
        return -1;
    }

    if (PAPI_OK != PAPI_query_event(PAPI_TOT_INS)) fprintf(stderr, "Cannot count PAPI_TOT_INS.\n");
    if (PAPI_OK != PAPI_query_event(PAPI_LST_INS)) fprintf(stderr, "Cannot count PAPI_LST_INS.\n");
    if (PAPI_OK != PAPI_query_event(PAPI_FP_INS )) fprintf(stderr, "Cannot count PAPI_FP_INS.\n" );
    if (PAPI_OK != PAPI_query_event(PAPI_L1_DCA )) fprintf(stderr, "Cannot count PAPI_L1_DCA.\n" );
    if (PAPI_OK != PAPI_query_event(PAPI_L1_DCM )) fprintf(stderr, "Cannot count PAPI_L1_DCM.\n" );
    if (PAPI_OK != PAPI_query_event(PAPI_L2_DCA )) fprintf(stderr, "Cannot count PAPI_L2_DCA.\n" );
    if (PAPI_OK != PAPI_query_event(PAPI_L2_DCM )) fprintf(stderr, "Cannot count PAPI_L2_DCM.\n" );
    if (PAPI_OK != PAPI_query_event(PAPI_TOT_CYC)) fprintf(stderr, "Cannot count PAPI_TOT_CYC.\n");
    long long countersOne[3];
    long long countersTwo[4];
    long long countersThree[1];

    for (int size = 0; size < 12; size++)
    {
        int half = sizes[size] / 2;

        double **A  = malloc (sizes[size] * sizeof(double *));
        double **B  = malloc (sizes[size] * sizeof(double *));
        double **C  = malloc (sizes[size] * sizeof(double *));
        double **D  = malloc (sizes[size] * sizeof(double *)); // used for checking validity.

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

        for (int i = 0; i < sizes[size]; i++)
        {
            A[i]  = malloc(sizes[size] * sizeof(double));
            B[i]  = malloc(sizes[size] * sizeof(double));
            C[i]  = malloc(sizes[size] * sizeof(double));
            D[i]  = malloc(sizes[size] * sizeof(double));

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

        for (int i = 0; i < sizes[size]; i++) 
        {
            for (int j = 0; j < sizes[size]; j++) 
            {
                A[i][j]  = rand_from(1, 3);
                B[i][j]  = rand_from(1, 3);
                C[i][j]  = 0;
            }
        }

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
            }
        }

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

        for (int i = 0; i < half; i++)
        {
            for (int j = 0; j < half; j++)
            {
                s1[i][j] = a21[i][j] + a22[i][j];
                s2[i][j] = s1[i][j]  - a11[i][j];
                s3[i][j] = a11[i][j] - a21[i][j];
                s4[i][j] = a12[i][j] -  s2[i][j];

                t1[i][j] = b12[i][j] - b11[i][j];
                t2[i][j] = b22[i][j] -  t1[i][j];
                t3[i][j] = b22[i][j] - b12[i][j];
                t4[i][j] = b21[i][j] -  t2[i][j];
            }
        }
        // printf("A is:\n");
        // for (int i = 0; i < sizes[size]; i++)
        // {
        //     for (int j = 0; j < sizes[size]; j++)
        //     {
        //         printf("%.2f ", A[i][j]);
        //     }
        //     printf("\n");
        // }

        // printf("B is:\n");
        // for (int i = 0; i < sizes[size]; i++)
        // {
        //     for (int j = 0; j < sizes[size]; j++)
        //     {
        //         printf("%.2f ", B[i][j]);
        //     }
        //     printf("\n");
        // }
        // struct timespec start_time, end_time;
        // long long elapsed_ns;
        // printf("\033[1;32mTiming data for winograd_strassen at size: %d\033[0m ", sizes[size]);
        // if (clock_gettime(CLOCK_REALTIME, &start_time) == -1) 
        // {
        //     printf("An error has occurred\n");
        //     return -1;
        // }
        printf("\033[1;32mTiming data for winograd_strassen at size: %d\033[0m", sizes[size]);
        if (PAPI_start_counters(PAPI_one, 3) != 0)
        {
            printf("AN ERROR HAS OCCURRED\n");
            return -1;
        }
        strassen(sizes[size], A, B, C, a11, a12, a21, a22, b11, b12, b21, b22, s1, s2, s3, s4, t1, t2, t3, t4, p1, p2, p3, p4, p5, p6, p7, u1, u2, u3, u4, u5, u6, u7);
        if (PAPI_stop_counters(countersOne, 3) != 0) 
        {
            printf("AN ERROR HAS OCCURRED\n");
            return -1;
        }
        print_one_data(countersOne, sizes[size]);
        reset_matrices(A, B, C, sizes[size]);

        if (PAPI_start_counters(PAPI_two, 4) != 0) 
        {
            printf("AN ERROR HAS OCCURRED\n");
            return -1;
        }
        strassen(sizes[size], A, B, C, a11, a12, a21, a22, b11, b12, b21, b22, s1, s2, s3, s4, t1, t2, t3, t4, p1, p2, p3, p4, p5, p6, p7, u1, u2, u3, u4, u5, u6, u7);
        if (PAPI_stop_counters(countersTwo, 4) != 0) 
        {
            printf("AN ERROR HAS OCCURRED\n");
            return -1;
        }
        print_two_data(countersTwo, sizes[size]);
        reset_matrices(A, B, C, sizes[size]);

        if (PAPI_start_counters(PAPI_three, 1) != 0) 
        {
            printf("AN ERROR HAS OCCURRED\n");
            return -1;
        }
        strassen(sizes[size], A, B, C, a11, a12, a21, a22, b11, b12, b21, b22, s1, s2, s3, s4, t1, t2, t3, t4, p1, p2, p3, p4, p5, p6, p7, u1, u2, u3, u4, u5, u6, u7);
        if (PAPI_stop_counters(countersThree, 1) != 0) 
        {
            printf("AN ERROR HAS OCCURRED\n");
            return -1;
        }
        print_three_data(countersThree, sizes[size]);

        // if (clock_gettime(CLOCK_REALTIME, &end_time) == -1) 
        // {
        //     printf("An error has occurred\n");
        //     return -1;
        // }
        // elapsed_ns = (end_time.tv_sec - start_time.tv_sec) * 1000000000LL + (end_time.tv_nsec - start_time.tv_nsec);
        // printf("CLOCK_REALTIME Clock elapsed: %lld ns\n", elapsed_ns);
        naive(A, B, D, sizes[size]);
        // for (int i = 0; i < sizes[size]; i++)
        // {
        //     for (int j = 0; j < sizes[size]; j++)
        //     {
        //         if (abs(C[i][j] - D[i][j]) != 0) 
        //         {
        //             printf("\033[1;31mError. C[%d][%d] is %f, \
        //             while D[%d][%d] is %.10f. \033[0m\n", i, j, C[i][j], i, j, D[i][j]);
        //         }
        //     }
        // }
    }

    return 0;
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
        // printf("a11 is:\n");
        // for (int i = 0; i < half; i++)
        // {
        //     for (int j = 0; j < half; j++)
        //     {
        //         printf("%.2f ", a11[i][j]);
        //     }
        //     printf("\n");
        // }
        // printf("\n");

        // printf("a12 is:\n");
        // for (int i = 0; i < half; i++)
        // {
        //     for (int j = 0; j < half; j++)
        //     {
        //         printf("%.2f ", a12[i][j]);
        //     }
        //     printf("\n");
        // }
        // printf("\n");

        // printf("a21 is:\n");
        // for (int i = 0; i < half; i++)
        // {
        //     for (int j = 0; j < half; j++)
        //     {
        //         printf("%.2f ", a21[i][j]);
        //     }
        //     printf("\n");
        // }
        // printf("\n");

        // printf("a22 is:\n");
        // for (int i = 0; i < half; i++)
        // {
        //     for (int j = 0; j < half; j++)
        //     {
        //         printf("%.2f ", a22[i][j]);
        //     }
        //     printf("\n");
        // }
        // printf("\n");

        // printf("b11 is:\n");
        // for (int i = 0; i < half; i++)
        // {
        //     for (int j = 0; j < half; j++)
        //     {
        //         printf("%.2f ", b11[i][j]);
        //     }
        //     printf("\n");
        // }
        // printf("\n");

        // printf("b12 is:\n");
        // for (int i = 0; i < half; i++)
        // {
        //     for (int j = 0; j < half; j++)
        //     {
        //         printf("%.2f ", b12[i][j]);
        //     }
        //     printf("\n");
        // }
        // printf("\n");

        // printf("b21 is:\n");
        // for (int i = 0; i < half; i++)
        // {
        //     for (int j = 0; j < half; j++)
        //     {
        //         printf("%.2f ", b21[i][j]);
        //     }
        //     printf("\n");
        // }
        // printf("\n");
        // printf("b22 is:\n");
        // for (int i = 0; i < half; i++)
        // {
        //     for (int j = 0; j < half; j++)
        //     {
        //         printf("%.2f ", b22[i][j]);
        //     }
        //     printf("\n");
        // }
        // calculate s1 --> s4 and t1 --> t4
        for (int i = 0; i < half; i++)
        {
            for (int j = 0; j < half; j++)
            {
                s1[i][j] = a21[i][j] + a22[i][j];
                s2[i][j] = s1[i][j]  - a11[i][j];
                s3[i][j] = a11[i][j] - a21[i][j];
                s4[i][j] = a12[i][j] -  s2[i][j];
                
                t1[i][j] = b12[i][j] - b11[i][j];
                t2[i][j] = b22[i][j] -  t1[i][j];
                t3[i][j] = b22[i][j] - b12[i][j];
                t4[i][j] = b21[i][j] -  t2[i][j];
            }
        }

        strassen(half, a11, b11,  p1, a11, a12, a21, a22, b11, b12, b21, b22, s1, s2, s3, s4, t1, t2, t3, t4, p1, p2, p3, p4, p5, p6, p7, u1, u2, u3, u4, u5, u6, u7); // p1 = a11 * b11
        strassen(half, a12, b21,  p2, a11, a12, a21, a22, b11, b12, b21, b22, s1, s2, s3, s4, t1, t2, t3, t4, p1, p2, p3, p4, p5, p6, p7, u1, u2, u3, u4, u5, u6, u7); // p2 = a12 * b21
        strassen(half, s1,   t1,  p3, a11, a12, a21, a22, b11, b12, b21, b22, s1, s2, s3, s4, t1, t2, t3, t4, p1, p2, p3, p4, p5, p6, p7, u1, u2, u3, u4, u5, u6, u7); // p3 = s1  * t1
        strassen(half, s2,   t2,  p4, a11, a12, a21, a22, b11, b12, b21, b22, s1, s2, s3, s4, t1, t2, t3, t4, p1, p2, p3, p4, p5, p6, p7, u1, u2, u3, u4, u5, u6, u7); // p4 = s2  * t2
        strassen(half, s3,   t3,  p5, a11, a12, a21, a22, b11, b12, b21, b22, s1, s2, s3, s4, t1, t2, t3, t4, p1, p2, p3, p4, p5, p6, p7, u1, u2, u3, u4, u5, u6, u7); // p5 = s3  * t3
        strassen(half, s4,   b22, p6, a11, a12, a21, a22, b11, b12, b21, b22, s1, s2, s3, s4, t1, t2, t3, t4, p1, p2, p3, p4, p5, p6, p7, u1, u2, u3, u4, u5, u6, u7); // p6 = s4  * b22
        strassen(half, a22,  t4,  p7, a11, a12, a21, a22, b11, b12, b21, b22, s1, s2, s3, s4, t1, t2, t3, t4, p1, p2, p3, p4, p5, p6, p7, u1, u2, u3, u4, u5, u6, u7); // p7 = a22 * t4
        
        // calculate u1 --> u7.
        for (int i = 0; i < half; i++) 
        {
            for (int j = 0; j < half; j++) 
            {
                u1[i][j] = p1[i][j] + p2[i][j];
                u2[i][j] = p1[i][j] + p4[i][j];
                u3[i][j] = u2[i][j] + p5[i][j];
                u4[i][j] = u3[i][j] + p7[i][j];
                u5[i][j] = u3[i][j] + p3[i][j];
                u6[i][j] = u2[i][j] + p3[i][j];
                u7[i][j] = u6[i][j] + p6[i][j];
            }
        }

        // set C11 --> C22.
        for (int i = 0; i < half; i++)
        {
            for (int j = 0; j < half; j++)
            {
                C[i][j]               = u1[i][j];  // C11
                C[i][j + half]        = u7[i][j];  // C12
                C[i + half][j]        = u4[i][j];  // C21
                C[i + half][j + half] = u5[i][j];  // C22
            }
        }

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
