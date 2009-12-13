#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "../lib/matrix.h"
/**********************************************************************/
int ID = -1;                /* Identify of MPI process */
int PROCESS_NUMBER = -1;    /* Number of MPI process */
#define CONSTANT_NUMBER 5
index_t tableOfConsts[CONSTANT_NUMBER];
#define MATRIX_SIZE               tableOfConsts[0]
#define MATRIX_ELEMENTS_NUMBER    tableOfConsts[1]
#define SUBMATRIX_SIZE            tableOfConsts[2]
#define SUBMATRIX_NUMBER          tableOfConsts[3]
#define SQRT_OF_SUBMATRIX_NUMBER  tableOfConsts[4]
/**********************************************************************/
typedef struct subMatrix {
	index_t    x;
	index_t    y;
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
	SQRT_OF_SUBMATRIX_NUMBER = --s;
	SUBMATRIX_NUMBER = SQRT_OF_SUBMATRIX_NUMBER * SQRT_OF_SUBMATRIX_NUMBER;
	MATRIX_SIZE = 0;
	MATRIX_ELEMENTS_NUMBER = 0;
	SUBMATRIX_SIZE = 0;
	return MPI_SUCCESS;
}

void placeSubMatrix(index_t x, index_t y, psubMatrix submatrix) {
	submatrix->x = x * SUBMATRIX_SIZE;
	submatrix->y = y * SUBMATRIX_SIZE;
	return;
}

inline void moveSubMatrixLeft(psubMatrix S) {
	S->x = (S->x + SUBMATRIX_SIZE) % MATRIX_SIZE;
	return;
}

inline void moveSubMatrixUp(psubMatrix S) {
	S->y =  (S->y + SUBMATRIX_SIZE) % MATRIX_SIZE;
	return;
}

inline unit_t *getSubMatrixElement(psubMatrix submatrix, index_t x, index_t y) {
	return getMatrixElement(submatrix->matrix, submatrix->x + x, submatrix->y + y);
}

inline unit_t getFromSubMatrix(const psubMatrix submatrix, index_t x, index_t y) {
	return *getSubMatrixElement(submatrix, x, y);
}

inline void setInSubMatrix(psubMatrix submatrix, index_t x, index_t y, unit_t value) {
	*getSubMatrixElement(submatrix, x, y) = value;
	return;
}

void subMatrixMul(const psubMatrix A, const psubMatrix B, psubMatrix AB) {
	int x, y, r;
	for(y = 0; y < SUBMATRIX_SIZE; y++)
		for(x = 0; x < SUBMATRIX_SIZE; x++)
			for(r = 0; r < SUBMATRIX_SIZE; r++)
				setInSubMatrix(AB, x, y, getFromSubMatrix(A, r, y) * getFromSubMatrix(B, x, r) + getFromSubMatrix(AB, x, y));
	return;
}

int sizeCheck(index_t rows, index_t columns) {
	if (rows != columns) {
		fprintf(stderr, "It is not a square matrix\n");
	   	return 1;
	}
	if (rows % SQRT_OF_SUBMATRIX_NUMBER != 0) {
	   	fprintf(stderr, "Size of the matrix read from file doesn't match process number\n");
		return 1;
	}
	if (MATRIX_SIZE == 0) {
		MATRIX_SIZE = rows;
		MATRIX_ELEMENTS_NUMBER = rows * rows;
		SUBMATRIX_SIZE = rows / SQRT_OF_SUBMATRIX_NUMBER;
	} else if (rows != MATRIX_SIZE) {
		fprintf(stderr, "Size of matrix read from file is doesn't match size of previously read matrix\n");
		return 1;
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
		close(file);
	} else {
		fprintf(stderr, "Can't open file '%s'\n", fileName);
		res = 1;
	}
	return res;
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

int Cannon_calculation(const pmatrix_t A, const pmatrix_t B, pmatrix_t AB, index_t max_shift) {
	subMatrix sA, sB, sAB;
	index_t x = ID % max_shift, y = ID / max_shift, shift;

	sA.matrix = A;
	sB.matrix = B;
	sAB.matrix = AB;
	placeSubMatrix((x+y) % max_shift, y, &sA);
	placeSubMatrix(x, (x+y) % max_shift, &sB);
	placeSubMatrix(x, y, &sAB);
	for(shift = 0; shift < max_shift; shift++) {
		subMatrixMul(&sA, &sB, &sAB);
		moveSubMatrixLeft(&sA);
		moveSubMatrixUp(&sB);
	}
	return MPI_Allreduce(AB->tab, A->tab, MATRIX_ELEMENTS_NUMBER, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}

#define MATRIX_A_FILE_NAME   argv[1]
#define MATRIX_B_FILE_NAME   argv[2]
#define MATRIX_AB_FILE_NAME  argv[3]

int main(int argc, char** argv) {
	pmatrix_t A, B, AB;

	if (mpiInit(&argc, &argv) != MPI_SUCCESS)
		fprintf(stderr, "Error in MPI initialization.\n"), exit(1);

	if ((ID == 0) && (argc < 4)) {
		printf("MPI Matrix multiplication (A * B = AB)\nUsage:\n\t%s <A matrix file (to read)> <B matrix file (to read)> <results matrix file (to write)>\n", argv[0]);
		MPI_Finalize(); return 1;
	}

	if (ID == 0) {
		if (readMatrixFromFile(MATRIX_A_FILE_NAME, &A) != 0)
			{ MPI_Finalize(); return 0; }
		if (readMatrixFromFile(MATRIX_B_FILE_NAME, &B) != 0)
			{ MPI_Finalize(); return 0; }
		printf("Data have been read. Starting calculation.\n");
	}

	if (MPI_Bcast(tableOfConsts, CONSTANT_NUMBER, MPI_INT, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
		fprintf(stderr, "[ID:%d] MPI_Bcast error.\n", ID), exit(1);

	if (initMatrix(&AB, MATRIX_SIZE, MATRIX_SIZE) == NULL)
		fprintf(stderr, "Can't allocate memory for AB matrix.\n"), exit(1);

	if (ID != 0) {
		if (initMatrix(&A, MATRIX_SIZE, MATRIX_SIZE) == NULL)
			fprintf(stderr, "Can't allocate memory for A matrix.\n"), exit(1);
		if (initMatrix(&B, MATRIX_SIZE, MATRIX_SIZE) == NULL)
			fprintf(stderr, "Can't allocate memory for B matrix.\n"), exit(1);
	}

	if (MPI_Bcast(A->tab, MATRIX_ELEMENTS_NUMBER, MPI_DOUBLE, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
		fprintf(stderr, "[ID:%d] MPI_Bcast error.\n", ID);
	if (MPI_Bcast(B->tab, MATRIX_ELEMENTS_NUMBER, MPI_DOUBLE, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
		fprintf(stderr, "[ID:%d] MPI_Bcast error.\n", ID);

	if (Cannon_calculation(A, B, AB, SQRT_OF_SUBMATRIX_NUMBER) != MPI_SUCCESS)
		fprintf(stderr, "[ID:%d] MPI_Allreduce error.\n", ID);

	if (ID == 0) {
		printf("Matrix multiplication has been completed.\n");
		writeFile(MATRIX_AB_FILE_NAME, A);
	}

	freeMatrix(&A);
	freeMatrix(&B);
	freeMatrix(&AB);

	MPI_Finalize();
	return 0;
}
/**********************************************************************/
