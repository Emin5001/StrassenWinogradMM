#ifndef MATRIX_H_
#define MATRIX_H_

#include <stdint.h>

#define MATRIX_ARITH_BEGIN __start_Matrix_Arith
#define MATRIX_ARITH_END   __stop_Matrix_Arith

#define MATRIX_INITIALIZER(X, ROW, COLUMN) \
        matrix_init(&(X) ,(ROW), (COLUMN))

#define autofree \
        __attribute__((cleanup(matrix_free)))

enum { NAIVE_ARITHMETIC, };

typedef struct _Matrix {
    uint32_t row;
    uint32_t column;
    double **values;
    
    double **(*New)(uint32_t row, uint32_t column);
    void  (*Delete)(void *);
} Matrix;

typedef struct _Matrix_Arith {
    Matrix (*Addition)(Matrix, const Matrix, const Matrix);
    Matrix (*Subtract)(Matrix, const Matrix, const Matrix);
    Matrix (*Multiply)(Matrix, const Matrix, const Matrix);
} Matrix_Arith;

void matrix_init(Matrix *, uint32_t, uint32_t);
void matrix_free(void *);

extern Matrix_Arith Naive_Matrix_Arith;

extern Matrix_Arith __start_Matrix_Arith[];
extern Matrix_Arith __stop_Matrix_Arith[];

#endif /* MATRIX_H_ */