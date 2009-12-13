#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "../lib/matrix.h"
/**********************************************************************/
int ID = 0;
int PROCESS_NUMBER = -1;
index_t tableOfConstants[] = { 0, 0, 0, 0, 0 } ;
#define TABLE_O_CONSTANS_NUMBER    4
#define MATRIX_A_ROWS_NUMBER       tableOfConstants[0] /* Number of matrix A rows */
#define MATRIX_A_COLUMNS_NUMBER    tableOfConstants[1] /* Number of matrix A columns and matrix B rows */
#define MATRIX_B_COLUMNS_NUMBER    tableOfConstants[2] /* Number of matrix B columns */
#define MESH_PROCESSOR_SIZE_NUMBER tableOfConstants[3] /* Size of processors mesh */
/**********************************************************************/
typedef struct subMatrix {
	index_t    x;
	index_t    y;
	index_t    rows;
	index_t    columns;
	pmatrix_t  matrix;
} subMatrix;

typedef subMatrix *psubMatrix;
/**********************************************************************/
int mpiInit(int *argc, char ***argv) {
	int s, s2, res;
	if ((res = MPI_Init(argc, argv)) != MPI_SUCCESS) return res;
	if ((res = MPI_Comm_rank(MPI_COMM_WORLD, &ID)) != MPI_SUCCESS) return res;
	if ((res = MPI_Comm_size(MPI_COMM_WORLD, &PROCESS_NUMBER)) != MPI_SUCCESS) return res;
	s2 = s = 1;
	while(s2 <= PROCESS_NUMBER) { s++; s2 = s * s; }
	MESH_PROCESSOR_SIZE_NUMBER = --s;
	return MPI_SUCCESS;
}

inline unit_t *getSubMatrixElement(psubMatrix A, index_t rows, index_t columns) {
	return getMatrixElement(A->matrix, A->y + rows, A->x + columns);
}

inline unit_t getFromSubMatrix(psubMatrix A, index_t rows, index_t columns) {
	return *getSubMatrixElement(A, rows, columns);
}

inline void setInSubMatrix(psubMatrix A, index_t rows, index_t columns, unit_t value) {
	*getSubMatrixElement(A, rows, columns) = value;
	return;
}

void subMatrixMul(const psubMatrix sA, const psubMatrix sB, psubMatrix sAB) {
	int x, y, r;

	for(y = 0; y < sA->rows; y++)
		for(x = 0; x < sB->columns; x++)
			for(r = 0; r < sA->columns; r++)
				setInSubMatrix(sAB, y, x, getFromSubMatrix(sA, y, r) * getFromSubMatrix(sB, r, x) + getFromSubMatrix(sAB, y, x));
	return;
}

int sizeCheck(index_t rows, index_t columns) {
	if (MATRIX_A_ROWS_NUMBER == 0) { /* Loading A matrix */
		MATRIX_A_ROWS_NUMBER = rows;
		MATRIX_A_COLUMNS_NUMBER = columns;
	} else { /* Loading B matrix */
		if (MATRIX_A_COLUMNS_NUMBER == rows) {
			MATRIX_B_COLUMNS_NUMBER = columns;
		}
		else {
			fprintf(stderr, "Matrix size don't match!\n Columns in B matrix should be equal to number of rows in A.\n");
			return 1;
		}
	}
	return 0;
}

int readMatrixFromFile(const char *fileName, pmatrix_t *matrix) {
	FILE *file; char res;
	printf("Reading matrix file: '%s'\n", fileName);
	if ((file = fopen(fileName, "r")) != NULL) {
		switch((res = readMatrix(file, matrix, &sizeCheck))) {
			case ERRORCODE_CANT_MALLOC:
				fprintf(stderr, "Can't allocate memory for matrix from file '%s'\n", fileName);
				break;
			case ERRORCODE_PARSER_ERROR:
				fprintf(stderr, "Syntax error in file '%s'\n", fileName);
				break;
			case ERRORCODE_SIZE_DONT_MATCH:
			case ERRORCODE_NOERRORS:
				break;
		}
		fclose(file);
	} else {
		fprintf(stderr, "Can't open file '%s'\n", fileName);
		res = 1;
	}
	return res;
}

index_t *lengthToStartIndex(const index_t *lengthsTab, const index_t size) {
	index_t i, *results = malloc(sizeof(index_t) * size);
	if (results == NULL) return NULL;
	results[0] = 0;
	for(i = 1; i < size; i++)
		results[i] = results[i -1] + lengthsTab[i -1];
	return results;
}

void printSubMatrix(psubMatrix s) {
	int x, y;
	for(x = 0; x < s->columns; x++) {
		for(y = 0; y < s->rows; y++)
			printf("%lf\t", getFromSubMatrix(s, x, y));
		printf("\n");
	}
}

#define ERROR_NO_ERRORS      0
#define ERROR_CANT_MALLOC    2
#define ERROR_MPI_ERROR      3

int SUMMA_calculation(const pmatrix_t A, const pmatrix_t B, pmatrix_t *AB, const index_t meshSideSize, const index_t aRowsLengths[], const index_t aColumnsLengths[], const index_t bColumnsLengths[]) {
	subMatrix sA, sB, sAB;
	matrix_t *tmpMatrix;
	index_t *aRowsBegins, *aColumnsBegins, *bColumnsBegins;
	index_t k, x = ID / meshSideSize, y = ID % meshSideSize;

	initMatrix(&tmpMatrix, (*AB)->rows, (*AB)->columns);
	aRowsBegins = lengthToStartIndex(aRowsLengths, meshSideSize);
	aColumnsBegins = lengthToStartIndex(aColumnsLengths, meshSideSize);
	bColumnsBegins = lengthToStartIndex(bColumnsLengths, meshSideSize);
	if (!(aRowsBegins && aColumnsBegins && bColumnsBegins))
		return ERROR_CANT_MALLOC;

	sA.matrix = A; sB.matrix = B; sAB.matrix = tmpMatrix;

	sAB.x = bColumnsBegins[x];
	sAB.y = aRowsBegins[y];
	sAB.rows = aRowsLengths[x];
	sAB.columns = bColumnsLengths[y];
	sA.y = bColumnsBegins[y];
	sA.rows = aRowsLengths[y];
	sB.x = bColumnsBegins[x];
	sB.columns = bColumnsLengths[x];

	for(k = 0; k < meshSideSize; k++) {
		sA.x = sB.y = aColumnsBegins[k];
		sA.columns = sB.rows = aColumnsLengths[k];
		subMatrixMul(&sA, &sB, &sAB);
	}

	free(aRowsBegins);
	free(aColumnsBegins);
	free(bColumnsBegins);
	if (MPI_Allreduce(tmpMatrix->tab, (*AB)->tab, (*AB)->rows * (*AB)->columns, MPI_UNIT, MPI_SUM, MPI_COMM_WORLD) != MPI_SUCCESS)
		return ERROR_MPI_ERROR;
	return ERROR_NO_ERRORS;
}

index_t *initLengthsTable(index_t size) {
	int i, pice, max;
	index_t *tmp = malloc(sizeof(index_t) * MESH_PROCESSOR_SIZE_NUMBER);
	if (tmp == NULL) return NULL;
	max = MESH_PROCESSOR_SIZE_NUMBER - 1;
	pice = size / MESH_PROCESSOR_SIZE_NUMBER;
	for(i = 0; i < max; i++)
		tmp[i] = pice;
	tmp[max] = size - pice * max;
	return tmp;
}

void writeFile(const char *fileName, pmatrix_t matrix) {
	FILE *file;
	file = fopen(fileName, "w");
	if (file == NULL) {
		fprintf(stderr, "Can't open file '%s' to write, matrix will be written on standard output\n", fileName);
		printMatrix(matrix);
	} else {
		printMatrixOnStream(matrix, file, UNIT_SPECIFIER " ");
		fclose(file);
		printf("Results matrix will be written to file '%s'.\n", fileName);
	}
	return;
}

#define PROGRAM_NAME             argv[0]
#define MATRIX_A_FILE_NAME       argv[1]
#define MATRIX_B_FILE_NAME       argv[2]
#define MATRIX_AB_FILE_NAME      argv[3]

int main(int argc, char** argv) {
	pmatrix_t A, B, AB;
	index_t *aRowsLengths, *aColumnsLengths, *bColumnsLengths;

	if (mpiInit(&argc, &argv) != MPI_SUCCESS)
		return 1;

	if (ID == 0) {
		if (argc < 4) {
			printf("MPI Scalable Universal Matrix Multiplication Algorithm implementation.\nUsage:\n\t%s <matrix A file (to read)> <matrix B file (to read)> <matrix result file (to write)>\n", PROGRAM_NAME);
			MPI_Abort(MPI_COMM_WORLD, 1); return 1;
		}
		if (readMatrixFromFile(MATRIX_A_FILE_NAME, &A) != 0)
			{ MPI_Abort(MPI_COMM_WORLD, 1); return 1; }
		if (readMatrixFromFile(MATRIX_B_FILE_NAME, &B) != 0)
			{ MPI_Abort(MPI_COMM_WORLD, 1); return 1; }
		printf("Data have been read. Now data will be send by broadcast.\n");
	}

	if (MPI_Bcast(tableOfConstants, TABLE_O_CONSTANS_NUMBER, MPI_INT, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
		fprintf(stderr, "[ID:%d] MPI_Bcast error.\n", ID);

	if (ID != 0) {
		initMatrix(&A, MATRIX_A_ROWS_NUMBER, MATRIX_A_COLUMNS_NUMBER);
		initMatrix(&B, MATRIX_A_COLUMNS_NUMBER, MATRIX_B_COLUMNS_NUMBER);
	}

	initMatrix(&AB, MATRIX_A_ROWS_NUMBER, MATRIX_B_COLUMNS_NUMBER);

	if (MPI_Bcast(A->tab, A->rows * A->columns, MPI_UNIT, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
		fprintf(stderr, "[ID:%d] MPI_Bcast error.\n", ID);

	if (MPI_Bcast(B->tab, B->rows * B->columns, MPI_UNIT, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
		fprintf(stderr, "[ID:%d] MPI_Bcast error.\n", ID);

	aRowsLengths = initLengthsTable(MATRIX_A_ROWS_NUMBER);
	aColumnsLengths = initLengthsTable(MATRIX_A_COLUMNS_NUMBER);
	bColumnsLengths = initLengthsTable(MATRIX_B_COLUMNS_NUMBER);

	if (!(aRowsLengths && aColumnsLengths && bColumnsLengths)) {
		fprintf(stderr, "Can't allocate lengths tables.\n");
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	SUMMA_calculation(A, B, &AB, MESH_PROCESSOR_SIZE_NUMBER, aRowsLengths, aColumnsLengths, bColumnsLengths);

	if (ID == 0) {
		printf("Matrix multiplication has been completed.\n");
		writeFile(MATRIX_AB_FILE_NAME, AB);
	}

	MPI_Finalize();
	return 0;
}
