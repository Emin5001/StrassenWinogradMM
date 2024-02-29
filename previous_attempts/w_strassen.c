#include <stdio.h>
#include <stdlib.h>
#include <time.h>
// #include <papi.h>
#include "matrix.h" 

int PAPI_measurement ();
void CLOCK_REALTIME_measurement ();
double rand_from      (double, double);
Matrix w_strassen     (Matrix, const Matrix, const Matrix, int);
void naive            (Matrix, const Matrix, const Matrix, int);
void print_matrix     (Matrix, char *);
void check_valid      (Matrix, Matrix, int);
void reset_matrices   (Matrix, Matrix, Matrix, int);
void print_one_data   (long long *, int);
void print_two_data   (long long *, int);
void print_three_data (long long *, int);

uint32_t sizes[12] = {2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096};

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
// }

// int PAPI_three[1] = { PAPI_TOT_CYC };

int main () 
{
    // PAPI_measurement();
    CLOCK_REALTIME_measurement();
    return 0;
}

void CLOCK_REALTIME_measurement ()
{
    int size = 8;
    Matrix matrixA, matrixB, matrixC, matrixD;

    srand(time (NULL));
    MATRIX_INITIALIZER (matrixA, size, size);
    MATRIX_INITIALIZER (matrixB, size, size);
    MATRIX_INITIALIZER (matrixC, size, size);
    MATRIX_INITIALIZER (matrixD, size, size);

    reset_matrices(matrixA, matrixB, matrixC, size);
    w_strassen(matrixC, matrixA, matrixB, size);
    print_matrix(matrixC, "Matrix C");
    naive (matrixD, matrixA, matrixB, size);
    print_matrix(matrixD, "Matrix D");
}

// int PAPI_measurement () 
// {
//     // if (PAPI_library_init(PAPI_VER_CURRENT) < 0) 
//     // {
//     //     fprintf(stderr, "Error intializing PAPI!");
//     //     return -1;
//     // }

//     // if (PAPI_OK != PAPI_query_event(PAPI_TOT_INS)) fprintf(stderr, "Cannot count PAPI_TOT_INS.\n");
//     // if (PAPI_OK != PAPI_query_event(PAPI_LST_INS)) fprintf(stderr, "Cannot count PAPI_LST_INS.\n");
//     // if (PAPI_OK != PAPI_query_event(PAPI_FP_INS )) fprintf(stderr, "Cannot count PAPI_FP_INS.\n" );
//     // if (PAPI_OK != PAPI_query_event(PAPI_L1_DCA )) fprintf(stderr, "Cannot count PAPI_L1_DCA.\n" );
//     // if (PAPI_OK != PAPI_query_event(PAPI_L1_DCM )) fprintf(stderr, "Cannot count PAPI_L1_DCM.\n" );
//     // if (PAPI_OK != PAPI_query_event(PAPI_L2_DCA )) fprintf(stderr, "Cannot count PAPI_L2_DCA.\n" );
//     // if (PAPI_OK != PAPI_query_event(PAPI_L2_DCM )) fprintf(stderr, "Cannot count PAPI_L2_DCM.\n" );
//     // if (PAPI_OK != PAPI_query_event(PAPI_TOT_CYC)) fprintf(stderr, "Cannot count PAPI_TOT_CYC.\n");
//     long long countersOne[3];
//     long long countersTwo[4];
//     long long countersThree[1];

//     Matrix matrixA, matrixB, matrixC, matrixD;

//     srand(time (NULL));
//     for (int size = 0; size < 12; size++) 
//     {
//         MATRIX_INITIALIZER (matrixA, sizes[size], sizes[size]);
//         MATRIX_INITIALIZER (matrixB, sizes[size], sizes[size]);
//         MATRIX_INITIALIZER (matrixC, sizes[size], sizes[size]);
//         MATRIX_INITIALIZER (matrixD, sizes[size], sizes[size]);

//         reset_matrices(matrixA, matrixB, matrixC, sizes[size]);

//         printf("\033[1;32mTiming data for winograd_strassen at size: %d\033[0m", sizes[size]);
//         if (PAPI_start_counters(PAPI_one, 3) != 0)
//         {
//             printf("AN ERROR HAS OCCURRED\n");
//             return -1;
//         }
//         w_strassen(matrixC, matrixA, matrixB, sizes[size]);
//         if (PAPI_stop_counters(countersOne, 3) != 0) 
//         {
//             printf("AN ERROR HAS OCCURRED\n");
//             return -1;
//         }
//         print_one_data(countersOne, sizes[size]);

//         reset_matrices(matrixA, matrixB, matrixC, sizes[size]);

//         if (PAPI_start_counters(PAPI_two, 4) != 0) 
//         {
//             printf("AN ERROR HAS OCCURRED\n");
//             return -1;
//         }
//         w_strassen(matrixC, matrixA, matrixB, sizes[size]);
//         if (PAPI_stop_counters(countersTwo, 4) != 0) 
//         {
//             printf("AN ERROR HAS OCCURRED\n");
//             return -1;
//         }
//         print_two_data(countersTwo, sizes[size]);

//         reset_matrices(matrixA, matrixB, matrixC, sizes[size]);

//         if (PAPI_start_counters(PAPI_three, 1) != 0) 
//         {
//             printf("AN ERROR HAS OCCURRED\n");
//             return -1;
//         }
//         w_strassen(matrixC, matrixA, matrixB, sizes[size]);
//         if (PAPI_stop_counters(countersThree, 1) != 0) 
//         {
//             printf("AN ERROR HAS OCCURRED\n");
//             return -1;
//         }
//         print_three_data(countersThree, sizes[size]);

//         matrixA.Delete(matrixA.values);
//         matrixB.Delete(matrixB.values);
//         matrixC.Delete(matrixC.values);
//     }

//     return -1;
// }

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

Matrix w_strassen(Matrix dest, const Matrix srcA, const Matrix srcB, int length) 
{
    Matrix_Arith *arith = &__start_Matrix_Arith[NAIVE_ARITHMETIC];
    
    if (length == 2) 
        return arith->Multiply(dest, srcA, srcB);

    int len = length / 2;

    autofree Matrix a11, a12, a21, a22, 
                    b11, b12, b21, b22,
                    c11, c12, c21, c22,
                    s1, s2, s3, s4,
                    t1, t2, t3, t4,
                    p1, p2 , p3 ,p4, p5, p6, p7,
                    u1, u2, u3, u4, u5, u6, u7;

    /* Initializer the matrix */
    MATRIX_INITIALIZER(a11, len, len); 
    MATRIX_INITIALIZER(a12, len, len); 
    MATRIX_INITIALIZER(a21, len, len); 
    MATRIX_INITIALIZER(a22, len, len);

    MATRIX_INITIALIZER(b11, len, len); 
    MATRIX_INITIALIZER(b12, len, len); 
    MATRIX_INITIALIZER(b21, len, len); 
    MATRIX_INITIALIZER(b22, len, len);

    MATRIX_INITIALIZER(c11, len, len); 
    MATRIX_INITIALIZER(c12, len, len); 
    MATRIX_INITIALIZER(c21, len, len); 
    MATRIX_INITIALIZER(c22, len, len);

    MATRIX_INITIALIZER(s1, len, len);
    MATRIX_INITIALIZER(s2, len, len);
    MATRIX_INITIALIZER(s3, len, len);
    MATRIX_INITIALIZER(s4, len, len);

    MATRIX_INITIALIZER(t1, len, len);
    MATRIX_INITIALIZER(t2, len, len);
    MATRIX_INITIALIZER(t3, len, len);
    MATRIX_INITIALIZER(t4, len, len);

    MATRIX_INITIALIZER(p1, len, len);  
    MATRIX_INITIALIZER(p2, len, len);  
    MATRIX_INITIALIZER(p3, len, len);  
    MATRIX_INITIALIZER(p4, len, len);
    MATRIX_INITIALIZER(p5, len, len);  
    MATRIX_INITIALIZER(p6, len, len); 
    MATRIX_INITIALIZER(p7, len, len);
    
    MATRIX_INITIALIZER(u1, len, len);  
    MATRIX_INITIALIZER(u2, len, len);  
    MATRIX_INITIALIZER(u3, len, len);  
    MATRIX_INITIALIZER(u4, len, len);
    MATRIX_INITIALIZER(u5, len, len);  
    MATRIX_INITIALIZER(u6, len, len); 
    MATRIX_INITIALIZER(u7, len, len);

    /* Divide matrix to four part */
    for (int i = 0; i < len; i++) {
        for (int j = 0; j < len; j++) {
            a11.values[i][j] = srcA.values[i][j];
            a12.values[i][j] = srcA.values[i][j + len];
            a21.values[i][j] = srcA.values[i + len][j];
            a22.values[i][j] = srcA.values[i + len][j + len];
            
            b11.values[i][j] = srcB.values[i][j];
            b12.values[i][j] = srcB.values[i][j + len];
            b21.values[i][j] = srcB.values[i + len][j];
            b22.values[i][j] = srcB.values[i + len][j + len];
        }
    }
    print_matrix(srcA, "A");
    print_matrix(a11, "a11");
    print_matrix(a12, "a12");
    print_matrix(a21, "a21");
    print_matrix(a22, "a22");
    /* Calculate seven formulas of Strassen Algorithm */
    arith->Addition(s1, a21, a22);
    arith->Subtract(s2, s1, a11 );
    arith->Subtract(s3, a11, a21);
    arith->Subtract(s4, a12, s2 );

    arith->Subtract(t1, b12, b11);
    arith->Subtract(t2, b22, t1 );
    arith->Subtract(t3, b22, b12);
    arith->Subtract(t4, b21, t2 );

    w_strassen(p1, a11, b11, len);
    w_strassen(p2, a12, b21, len);
    w_strassen(p3, s1, t1, len  );
    w_strassen(p4, s2, t2, len  );
    w_strassen(p5, s3, t3, len  );
    w_strassen(p6, s4, b22, len );
    w_strassen(p7, a22, t4, len );

    arith->Addition(u1, p1, p2  );
    arith->Addition(u2, p1, p4  );
    arith->Addition(u3, u2, p5  );
    arith->Addition(u4, u3, p7  );
    arith->Addition(u5, u3, p3  );
    arith->Addition(u6, u2, p3  );
    arith->Addition(u7, u6, p6  );

    arith->Addition(c11, p1, p2 );
    arith->Addition(c21, u3, p7 );
    arith->Addition(c22, u3, p3 );
    arith->Addition(c12, u6, p6 );
    
    /* Store the answer of matrix multiplication */
    for (int i = 0; i < len; i++) {
        for (int j = 0; j < len; j++) {
            dest.values[i][j]              = c11.values[i][j];
            dest.values[i][j + len]        = c12.values[i][j];
            dest.values[i + len][j]        = c21.values[i][j];
            dest.values[i + len][j + len]  = c22.values[i][j];
        }
    }

    return dest;
}

void naive(Matrix res, Matrix a, Matrix b, int size) 
{
    double sum = 0.0;
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            sum = 0.0;
            for (int k = 0; k < size; k++) {
                sum += a.values[i][k] * b.values[k][j];
            }
            res.values[i][j] += sum;
        }
    }
}

double rand_from (double min, double max) 
{
    double rand_decimal = ((double) rand() / RAND_MAX) * (max - min) + min;
    
    // rounds to nearest 2 decimal points.
    return (int)(rand_decimal * 100.0 + 0.5) / 100.0;
}

void print_matrix (Matrix m, char *identifier) 
{
    printf("Now printing out %s\n", identifier);

    for (int r = 0; r < sizes[0]; r++) {
        for (int c = 0; c < sizes[0]; c++) {
            printf("%.2f ", m.values[r][c]);
        }
        printf("\n");
    }
}

void check_valid (Matrix a, Matrix b, int len) {
    for (int i = 0; i < len; i++) {
        for (int j = 0; j < len; j++) {
            if (abs(a.values[i][j] - b.values[i][j]) != 0.0) {
                printf("\033[1;31mError. a.values[%d][%d] is %f, \
                while b.values[%d][%d] is %.10f. \033[0m\n", i, j,  \
                a.values[i][j], i, j, b.values[i][j]);
            }
        }
    }
}

void reset_matrices(Matrix a, Matrix b, Matrix c, int size) 
{
    for (int i = 0; i < size; i++) 
    {
        for (int j = 0; j < size; j++) 
        {
            a.values[i][j] = i + j;
            b.values[i][j] = i + j;
            c.values[i][j] = 0.0;
        }
    }
}
