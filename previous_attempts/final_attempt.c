#include <stdio.h>
#include <stdlib.h>
#include <papi.h>
#include <time.h>
#include <math.h>

typedef struct {
    double **p1, **p2, **p3, **p4, **p5, **p6, **p7;
    double **c11, **c12, **c21, **c22;
    double **s1, **s2, **s3, **s4, **t1, **t2, **t3, **t4;
    double **a11, **a12, **a21, **a22, **b11, **b12, **b21, **b22;
} StrassenAuxMatrices;

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


void strassen              (int, double **, double **, double **, int, StrassenAuxMatrices *);
void naive                 (double **, double **, double **, int);
void matrix_add            (double **, double **, double **, int);
void matrix_subtr          (double **, double **, double **, int);
void reset_matrices        (double **, double **, double **, int);
void initializeAuxMatrices (int, StrassenAuxMatrices *, int);
void print_one_data        (long long *, int);
void print_two_data        (long long *, int);
void print_three_data      (long long *, int);
double rand_from           (double, double);
void print_matrix          (double **, int);
double** allocateMatrix    (int);


int main()
{
    if (PAPI_library_init (PAPI_VER_CURRENT) < 0)
    {
        fprintf(stderr, "Error initializing PAPI!");
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

    int sizes[] = {2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048};
    for (int s = 0; s < 11; s++)
    {
        int size = sizes[s];
        int log_n = (int) (log (size) / log (2));

        StrassenAuxMatrices aux[log_n];
        initializeAuxMatrices(size, aux, log_n);

        double **A = (double **) malloc (size * sizeof (double *));
        double **B = (double **) malloc (size * sizeof (double *));
        double **C = (double **) malloc (size * sizeof (double *));
        double **D = (double **) malloc (size * sizeof (double *));

        for (int i = 0; i < size; i++)
        {
            A[i] = (double *) malloc (size * sizeof (double));
            B[i] = (double *) malloc (size * sizeof (double));
            C[i] = (double *) malloc (size * sizeof (double));
            D[i] = (double *) malloc (size * sizeof (double));
        }

        for (int i = 0; i < size; i++)
        {
            for (int j = 0; j < size; j++)
            {
                A[i][j] = rand_from (1, 3);
                B[i][j] = rand_from (1, 3);
            }
        }


        struct timespec start_time, end_time;
        long long elapsed_ns;
        printf("\033[1;32mTiming data for winograd_strassen at size: %d\033[0m ", size);
        if (clock_gettime(CLOCK_REALTIME, &start_time) == -1) 
        {
            printf("An error has occurred\n");
            return -1;
        }
        strassen (size, A, B, C, 0, aux);
        if (clock_gettime(CLOCK_REALTIME, &end_time) == -1) 
        {
            printf("An error has occurred\n");
            return -1;
        }

        elapsed_ns = (end_time.tv_sec - start_time.tv_sec) * 1000000000LL + (end_time.tv_nsec - start_time.tv_nsec);
        printf("CLOCK_REALTIME Clock elapsed: %lld ns\n", elapsed_ns);

        // printf("\033[1;32mTiming data for winograd_strassen at size: %d\033[0m", size);
        // if (PAPI_start_counters(PAPI_one, 3) != 0)
        // {
        //     printf("AN ERROR HAS OCCURRED\n");
        //     return -1;
        // }
        // strassen(size, A, B, C, 0, aux);
        // naive (A, B, D, size);
        // if (PAPI_stop_counters(countersOne, 3) != 0) 
        // {
        //     printf("AN ERROR HAS OCCURRED\n");
        //     return -1;
        // }
        // print_one_data(countersOne, size);
        // reset_matrices(A, B, C, size);

        // if (PAPI_start_counters(PAPI_two, 4) != 0) 
        // {
        //     printf("AN ERROR HAS OCCURRED\n");
        //     return -1;
        // }
        // strassen(size, A, B, C, 0, aux);
        // if (PAPI_stop_counters(countersTwo, 4) != 0) 
        // {
        //     printf("AN ERROR HAS OCCURRED\n");
        //     return -1;
        // }
        // print_two_data(countersTwo, size);
        // reset_matrices(A, B, C, size);

        // if (PAPI_start_counters(PAPI_three, 1) != 0) 
        // {
        //     printf("AN ERROR HAS OCCURRED\n");
        //     return -1;
        // }
        // strassen(size, A, B, C, 0, aux);
        // if (PAPI_stop_counters(countersThree, 1) != 0) 
        // {
        //     printf("AN ERROR HAS OCCURRED\n");
        //     return -1;
        // }
        // print_three_data(countersThree, size);

        // naive (A, B, D, size);
    }

    return 0;
}

void strassen(int n, double **A, double **B, double **C, int depth, StrassenAuxMatrices *aux)
{
    if (n == 2) {
        C[0][0] = (A[0][0] * B[0][0]) + (A[0][1] * B[1][0]);
        C[0][1] = (A[0][0] * B[0][1]) + (A[0][1] * B[1][1]);
        C[1][0] = (A[1][0] * B[0][0]) + (A[1][1] * B[1][0]);
        C[1][1] = (A[1][0] * B[0][1]) + (A[1][1] * B[1][1]);
        return;
    }

    int newSize = n / 2;

    for (int i = 0; i < newSize; i++)
    {
        for (int j = 0; j < newSize; j++)
        {
            aux[depth].a11[i][j] = A[i][j];
            aux[depth].a12[i][j] = A[i][j + newSize];
            aux[depth].a21[i][j] = A[i + newSize][j];
            aux[depth].a22[i][j] = A[i + newSize][j + newSize];

            aux[depth].b11[i][j] = B[i][j];
            aux[depth].b12[i][j] = B[i][j + newSize];
            aux[depth].b21[i][j] = B[i + newSize][j];
            aux[depth].b22[i][j] = B[i + newSize][j + newSize];
        }
    }

    matrix_add   (aux[depth].a12, aux[depth].a22, aux[depth].s1, newSize);
    matrix_subtr (aux[depth].s1,  aux[depth].a11, aux[depth].s2, newSize);
    matrix_subtr (aux[depth].a11, aux[depth].a21, aux[depth].s3, newSize);
    matrix_subtr (aux[depth].a12, aux[depth].s2 , aux[depth].s4, newSize);

    matrix_subtr (aux[depth].b12, aux[depth].b11, aux[depth].t1, newSize);
    matrix_subtr (aux[depth].b22, aux[depth].t1 , aux[depth].t2, newSize);
    matrix_subtr (aux[depth].b22, aux[depth].b12, aux[depth].t3, newSize);
    matrix_subtr (aux[depth].b21, aux[depth].t2 , aux[depth].t4, newSize);

    strassen (newSize, aux[depth].a11, aux[depth].b11, aux[depth].p1, depth + 1, aux);
    strassen (newSize, aux[depth].a12, aux[depth].b21, aux[depth].p2, depth + 1, aux);
    strassen (newSize, aux[depth].s1,  aux[depth].t1,  aux[depth].p3, depth + 1, aux);
    strassen (newSize, aux[depth].s2,  aux[depth].t2,  aux[depth].p4, depth + 1, aux);
    strassen (newSize, aux[depth].s3,  aux[depth].t3,  aux[depth].p5, depth + 1, aux);
    strassen (newSize, aux[depth].s4,  aux[depth].b22, aux[depth].p6, depth + 1, aux);
    strassen (newSize, aux[depth].a22, aux[depth].t4 , aux[depth].p7, depth + 1, aux);


    for (int i = 0; i < newSize; i++)
    {
        for (int j = 0; j < newSize; j++)
        {
            aux[depth].c11[i][j] = aux[depth].p1[i][j] + aux[depth].p2[i][j];
            aux[depth].c12[i][j] = aux[depth].p1[i][j] + aux[depth].p4[i][j] + aux[depth].p3[i][j] + aux[depth].p6[i][j]; 
            aux[depth].c21[i][j] = aux[depth].p1[i][j] + aux[depth].p4[i][j] + aux[depth].p5[i][j] + aux[depth].p7[i][j]; 
            aux[depth].c22[i][j] = aux[depth].p1[i][j] + aux[depth].p4[i][j] + aux[depth].p5[i][j] + aux[depth].p3[i][j];  
        }
    }

    for (int i = 0; i < newSize; i++)
    {
        for (int j = 0; j < newSize; j++)
        {
            C[i][j]                     = aux[depth].c11[i][j];
            C[i][j + newSize]           = aux[depth].c12[i][j];
            C[i + newSize][j]           = aux[depth].c21[i][j];
            C[i + newSize][j + newSize] = aux[depth].c22[i][j];
        }
    }
}

double** allocateMatrix (int size)
{
    double **matrix = (double **) malloc (size * sizeof (double *));

    for (int i = 0; i < size; i++)
    {
        matrix[i] = (double *) malloc (size * sizeof (double));
    }

    return matrix;
}

void initializeAuxMatrices(int maxN, StrassenAuxMatrices *aux, int LOG_N) {
    int size = maxN;
    for (int depth = 0; depth < LOG_N; depth++) {
        aux[depth].a11 = allocateMatrix(size);
        aux[depth].a12 = allocateMatrix(size);
        aux[depth].a21 = allocateMatrix(size);
        aux[depth].a22 = allocateMatrix(size);

        aux[depth].b11 = allocateMatrix(size);
        aux[depth].b12 = allocateMatrix(size);
        aux[depth].b21 = allocateMatrix(size);
        aux[depth].b22 = allocateMatrix(size);

        aux[depth].c11 = allocateMatrix(size);
        aux[depth].c12 = allocateMatrix(size);
        aux[depth].c21 = allocateMatrix(size);
        aux[depth].c22 = allocateMatrix(size);

        aux[depth].s1 = allocateMatrix(size);
        aux[depth].s2 = allocateMatrix(size);
        aux[depth].s3 = allocateMatrix(size);
        aux[depth].s4 = allocateMatrix(size);

        aux[depth].t1 = allocateMatrix(size);
        aux[depth].t2 = allocateMatrix(size);
        aux[depth].t3 = allocateMatrix(size);
        aux[depth].t4 = allocateMatrix(size);

        aux[depth].p1 = allocateMatrix(size);
        aux[depth].p2 = allocateMatrix(size);
        aux[depth].p3 = allocateMatrix(size);
        aux[depth].p4 = allocateMatrix(size);
        aux[depth].p5 = allocateMatrix(size);
        aux[depth].p6 = allocateMatrix(size);
        aux[depth].p7 = allocateMatrix(size);

        size /= 2;  // Halve the size for the next depth
    }
}

void matrix_add (double **a, double **b, double **res, int size)
{
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            res[i][j] = a[i][j] + b[i][j];
        }
    }
}

void matrix_subtr (double **a, double **b, double **res, int size)
{
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            res[i][j] = a[i][j] - b[i][j];
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