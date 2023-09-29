#include <stdlib.h>
#include "matrix.h"

static double **matrix_alloc(uint32_t row, uint32_t column)
{
    double **matrix = (double **) calloc(row + row * column, sizeof(**matrix));
    double *matrixTemp = (double *) (matrix + row);

    for (uint32_t i = 0; i < row; i++) {
        matrix[i] = matrixTemp;
        matrixTemp += column;
    }

    return matrix;
}

__attribute__((always_inline)) 
inline void matrix_free(void *ptr)
{
    Matrix *matrix = (Matrix *) ptr;
    return free(matrix->values);
}

static Matrix matrix_add(Matrix res, const Matrix a, const Matrix b)
{
    for (uint32_t i = 0; i < res.row; i++)
        for (uint32_t j = 0; j < res.column; j++)
            res.values[i][j] = a.values[i][j] + b.values[i][j];

    return res;
}

static Matrix matrix_subtract(Matrix res, const Matrix a, const Matrix b)
{
    for (uint32_t i = 0; i < res.row; i++)
        for (uint32_t j = 0; j < res.column; j++)
            res.values[i][j] = a.values[i][j] - b.values[i][j];

    return res;
}

static Matrix matrix_multiply(Matrix res, const Matrix a, const Matrix b)
{
    // TODO: CHANGE THIS.
    // double S1 = a.values[1][0] + a.values[1][1];
    // double S2 = S1 - a.values[0][0];
    // double S3 = a.values[0][0] - a.values[1][0];
    // double S4 = a.values[0][1] - S2;

    // double T1 = b.values[0][1] - b.values[0][0];
    // double T2 = b.values[1][1] - T1;
    // double T3 = b.values[1][1] - b.values[0][1];
    // double T4 = b.values[1][0] - T2;

    // double P1 = a.values[0][0] * b.values[0][0];
    // double P2 = a.values[0][1] * b.values[1][0];
    // double P3 = S1 * T1;
    // double P4 = S2 * T2;
    // double P5 = S3 * T3;
    // double P6 = S4 * b.values[1][1];
    // double P7 = a.values[1][1] * T4;

    // double U1 = P1 + P1;
    // double U2 = P1 + P4;
    // double U3 = U2 + P5;
    // double U4 = U3 + P7;
    // double U5 = U3 + P3;
    // double U6 = U2 + P3;
    // double U7 = U6 + P6;

    // res.values[0][0] = U1;
    // res.values[0][1] = U7;
    // res.values[1][0] = U4;
    // res.values[1][1] = U5;
  
    // return res;

    double m1 = (a.values[0][0] + a.values[1][1]) * (b.values[0][0] + b.values[1][1]);
    double m2 = (a.values[1][0] + a.values[1][1]) *  b.values[0][0];
    double m3 = (b.values[0][1] - b.values[1][1]) *  a.values[0][0];
    double m4 = (b.values[1][0] - b.values[0][0]) *  a.values[1][1];
    double m5 = (a.values[0][0] + a.values[0][1]) *  b.values[1][1];
    double m6 = (a.values[1][0] - a.values[0][0]) * (b.values[0][0] + b.values[0][1]);
    double m7 = (a.values[0][1] - a.values[1][1]) * (b.values[1][0] + b.values[1][1]);

    res.values[0][0] = m1 + m4 - m5 + m7;
    res.values[0][1] = m3 + m5;
    res.values[1][0] = m2 + m4;
    res.values[1][1] = m1 - m2 + m3 + m6;
  
    return res;
}

void matrix_init(Matrix *matrix, uint32_t row, uint32_t column)
{
    matrix->row    = row;
    matrix->column = column;
    matrix->values = NULL;
    matrix->New    = matrix_alloc;
    matrix->Delete = free;

    matrix->values = matrix->New(row, column);
}

Matrix_Arith Naive_Matrix_Arith __attribute__((section("Matrix_Arith"))) = {
    .Addition = matrix_add,
    .Subtract = matrix_subtract,
    .Multiply = matrix_multiply,
};