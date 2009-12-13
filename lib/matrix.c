#include "matrix.h"

pmatrix_t initMatrix(pmatrix_t *matrix, index_t rows, index_t columns) {
	index_t y;
	*matrix = malloc(sizeof(matrix_t));
	if (matrix == NULL) return NULL;
	(*matrix)->rows = rows;
	(*matrix)->columns = columns;
	(*matrix)->tab = malloc(sizeof(unit_t) * rows * columns);
	if ((*matrix)->tab == NULL) return NULL;
	(*matrix)->ptab = malloc(sizeof(unit_t *) * rows);
	if ((*matrix)->ptab == NULL) return NULL;
	for(y = 0; y < rows; y++)
		(*matrix)->ptab[y] = (*matrix)->tab + y * columns;
	return (*matrix);
}

inline pmatrix_t initCopyOfMatrix(pmatrix_t *matrix, const pmatrix_t original) {
	return initMatrix(matrix, original->rows, original->columns);
}

void freeMatrix(pmatrix_t *matrix) {
	free((*matrix)->ptab);
	free((*matrix)->tab);
	free(*matrix);
}

inline unit_t *getMatrixElement(pmatrix_t matrix, index_t rows, index_t columns) {
	return &(matrix->ptab[rows][columns]);
}

inline unit_t getFromMatrix(pmatrix_t matrix, index_t rows, index_t columns) {
	return *getMatrixElement(matrix, rows, columns);
}

inline void setInMatrix(pmatrix_t matrix, index_t rows, index_t columns, unit_t value) {
	*getMatrixElement(matrix, rows, columns) = value;
	return;
}

inline index_t getMatrixElementNumber(pmatrix_t matrix) {
	return (matrix->rows * matrix->columns);
}

void printMatrixOnStream(pmatrix_t matrix, FILE * stream, const char *unit_t_format) {
	index_t x, y;
	fprintf(stream, "%d %d\n", matrix->rows, matrix->columns);
	for(y = 0; y < matrix->rows; y++) {
		for(x = 0; x < matrix->columns; x++)
			fprintf(stream, unit_t_format, getFromMatrix(matrix, y, x));
		fprintf(stream, "\n");
	}
}

inline void printMatrix(pmatrix_t matrix) {
	printMatrixOnStream(matrix, stdout, UNIT_SPECIFIER "\t");
}

static int readFile(FILE *file, pmatrix_t *matrix, index_t rows, index_t columns) {
	index_t x, y;
	initMatrix(matrix, rows, columns);
	if (matrix == NULL)
		return ERRORCODE_CANT_MALLOC;
	for(y = 0; y < rows; y++) {
		for(x = 0; x < columns; x++)
			if (fscanf(file, UNIT_SPECIFIER, getMatrixElement(*matrix, y, x)) != 1)
				return ERRORCODE_PARSER_ERROR; 
	}
	return ERRORCODE_NOERRORS;
}

int readMatrix(FILE *file, pmatrix_t *matrix, pcheck_size_function checkSize) {
	index_t rows, columns, x, y;
	if (fscanf(file, "%d %d\n", &rows, &columns) == 2) {
		if (checkSize(rows, columns))
			return ERRORCODE_SIZE_DONT_MATCH;
		return readFile(file, matrix, rows, columns);
	} else
		return ERRORCODE_PARSER_ERROR;
}


