#ifndef __MATRIX_H__
#define __MATRIX_H__
#include <stdio.h>
#include <stdlib.h>
#include "unit.h"
#include "error.h"

typedef struct matrix_t {
	index_t    rows;
	index_t    columns;
	unit_t    *tab;
	unit_t   **ptab;
} matrix_t;

typedef matrix_t *pmatrix_t;
typedef int (*pcheck_size_function)(index_t, index_t);

pmatrix_t          initMatrix(pmatrix_t *matrix, index_t rows, index_t columns);
inline pmatrix_t   initCopyOfMatrix(pmatrix_t *matrix, const pmatrix_t original);
void               freeMatrix(pmatrix_t *matrix);
inline unit_t      *getMatrixElement(pmatrix_t matrix, index_t rows, index_t columns);
inline unit_t      getFromMatrix(pmatrix_t matrix, index_t rows, index_t columns);
inline void        setInMatrix(pmatrix_t matrix, index_t rows, index_t columns, unit_t value);
inline index_t     getMatrixElementNumber(pmatrix_t matrix);
void               printMatrixOnStream(pmatrix_t matrix, FILE * stream, const char *unit_t_format);
inline void        printMatrix(pmatrix_t matrix);
int                readMatrix(FILE *file, pmatrix_t *matrix, pcheck_size_function checkSize);

#endif
