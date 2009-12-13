#include <stdio.h>
#include <string.h>
#include <mpi.h>
#include "../lib/matrix.h"
#include "../lib/vector.h"

#define PARSER_BUFFER 20

#define INFINITY         9999999
#define INFINITY_IN_FILE "-"

int ID = -1;
int PROCESS_NUMBER = -1;
int NODES_NUMBER = -1;

index_t calculateLength(index_t id) {
	if (id < (PROCESS_NUMBER -1))
		return NODES_NUMBER / PROCESS_NUMBER;
	else
		return NODES_NUMBER - (id * (NODES_NUMBER / PROCESS_NUMBER));
}

int readGraphLine(FILE *file, pmatrix_t lenMatrix, index_t y) {
	char tmpRead[PARSER_BUFFER];
	index_t x;

	for(x = 0; x < lenMatrix->rows; x++) {
		if (fscanf(file, "%s", tmpRead) != 1)
			return ERRORCODE_PARSER_ERROR;
		if (strcmp(tmpRead, INFINITY_IN_FILE) != 0)
			setInMatrix(lenMatrix, y, x, atof(tmpRead));
		else
			setInMatrix(lenMatrix, y, x, INFINITY);
	}
	return ERRORCODE_NOERRORS;
}

int readGraphFile(char *fileName, pmatrix_t *lenMatrix) {
	FILE *file;
	index_t res, y;

	printf("Reading graph file: '%s'\n", fileName);
	if ((file = fopen(fileName, "r")) != NULL) {
		if (fscanf(file, "%d", &NODES_NUMBER) == 1) {
			if (initMatrix(lenMatrix, NODES_NUMBER, NODES_NUMBER) == NULL)
				return ERRORCODE_CANT_MALLOC;

			for(y = 0; y < NODES_NUMBER; y++)
				if ((res = readGraphLine(file, *lenMatrix, y)) != ERRORCODE_NOERRORS)
					return res;
		} else
			return ERRORCODE_PARSER_ERROR;
		fclose(file);
	} else
		return ERRORCODE_CANT_FILE_OPEN;
	return ERRORCODE_NOERRORS;
}

void writeGraphFile(FILE *file, pmatrix_t lenMatrix) {
	index_t x, y;
	unit_t tmp;

	fprintf(file, "%d\n", NODES_NUMBER);
	for(y = 0; y < NODES_NUMBER; y++) {
		for(x = 0; x < NODES_NUMBER; x++) {
			tmp = getFromMatrix(lenMatrix, y, x);
			if (tmp < INFINITY)
				fprintf(file, UNIT_SPECIFIER "\t", tmp);
			else
				fprintf(file, "-\t");
		}
		fprintf(file, "\n");
	}
}

void writeFile(char *fileName, pmatrix_t lenMatrix) {
	FILE *file;
	if ((file = fopen(fileName, "w")) != NULL) {
		writeGraphFile(file, lenMatrix);
		fclose(file);
	} else {
		fprintf(stderr, "Can't open file '%s' to write, matrix will be written on standard output\n", fileName);
		writeGraphFile(stdout, lenMatrix);
	}
}

int Dijkstra_calculation(pmatrix_t lenMatrix, index_t vertex, unit_t *d) {
	index_t i, k, mini;
	pvector_t p;

	if (initVector(&p, NODES_NUMBER) == NULL)
		return ERRORCODE_CANT_MALLOC;

	for (i = 0; i < NODES_NUMBER; i++) {
		d[i] = INFINITY;
		setInVector(p, i, 0); /* the i-th element has not yet been visited */
	}
	d[vertex] = 0;

	for(k = 0; k < NODES_NUMBER; k++) {
		mini = -1;
		for (i = 0; i < NODES_NUMBER; i++)
			if (!getFromVector(p, i) && ((mini == -1) || (d[i] < d[mini])))
				mini = i;
			setInVector(p, mini, 1);
			for (i = 0; i < NODES_NUMBER; ++i)
				if (getFromMatrix(lenMatrix, mini, i))
					if (d[mini] + getFromMatrix(lenMatrix, mini, i) < d[i])
						d[i] = d[mini] + getFromMatrix(lenMatrix, mini, i);
	}
	return ERRORCODE_NOERRORS;
}

int mpiInit(int *argc, char ***argv) {
	int res;
	if ((res = MPI_Init(argc, argv)) != MPI_SUCCESS) return res;
	if ((res = MPI_Comm_rank(MPI_COMM_WORLD, &ID)) != MPI_SUCCESS) return res;
	if ((res = MPI_Comm_size(MPI_COMM_WORLD, &PROCESS_NUMBER)) != MPI_SUCCESS) return res;
	return res;
}

int AbortAndExit(int result, char *reason) {
	fprintf(stderr, "%s\n", reason);
	MPI_Abort(MPI_COMM_WORLD, result);
	return result;
}

#define ARGV_PROGRAM_NAME argv[0]
#define ARGV_FILE_NAME    argv[1]
#define ARGV_RESULT_NAME  argv[2]

int main(int argc, char** argv) {
	pmatrix_t lenMatrix, result;
	index_t i, startVertex;
	MPI_Status status;

	if (argc < 3) {
		fprintf(stderr, "MPI Dijstra Parallel Algorithm implementation.\nUsage: \n\t%s <graph file (to read)> <result file (to write)>\n", argv[0]);
		return 1;
	}

	if (mpiInit(&argc, &argv) != MPI_SUCCESS)
		return AbortAndExit(ERRORCODE_MPI_ERROR, "Can't invoke MPI");

	if ((ID == 0) && ((i = readGraphFile(ARGV_FILE_NAME, &lenMatrix)) != ERRORCODE_NOERRORS))
		return AbortAndExit(ERRORCODE_MPI_ERROR, "Error while reading file");

	if (MPI_Bcast(&NODES_NUMBER, 1, MPI_INT, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
		return AbortAndExit(ERRORCODE_MPI_ERROR, "MPI_Bcast error");

	if ((ID != 0) && (initMatrix(&lenMatrix, NODES_NUMBER, NODES_NUMBER) == NULL))
		return AbortAndExit(ERRORCODE_CANT_MALLOC, "Cannot allocate memory");

	if (MPI_Bcast(lenMatrix->tab, getMatrixElementNumber(lenMatrix), MPI_UNIT, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
		return AbortAndExit(ERRORCODE_MPI_ERROR, "MPI_Bcast error");

	if (initMatrix(&result, (ID == 0) ? NODES_NUMBER : calculateLength(ID), NODES_NUMBER) == NULL)
		return AbortAndExit(ERRORCODE_CANT_MALLOC, "Cannot allocate memory");

	if (ID == 0) printf("Calculation start\n");
	startVertex = ID * (NODES_NUMBER/PROCESS_NUMBER);
	for(i = 0; i < calculateLength(ID); i++) {
		if (Dijkstra_calculation(lenMatrix, startVertex + i, result->ptab[i]) != ERRORCODE_NOERRORS)
			return AbortAndExit(ERRORCODE_MPI_ERROR, "Error in calculation");
	}

	if (ID == 0) {
		startVertex = calculateLength(0);
		for(i = 1; i < PROCESS_NUMBER; i++) {
			if (MPI_Recv(result->ptab[startVertex], calculateLength(i) * NODES_NUMBER, MPI_UNIT, i, 0, MPI_COMM_WORLD, &status) != MPI_SUCCESS)
				return AbortAndExit(ERRORCODE_MPI_ERROR, "MPI_Recv error");
			startVertex += calculateLength(i);
		}
		printf("Calculation finished. Results will be written to file: '%s'\n", ARGV_RESULT_NAME);
		writeFile(ARGV_RESULT_NAME, result);
	} else
		if (MPI_Send(result->tab, NODES_NUMBER * calculateLength(ID), MPI_UNIT, 0, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
			return AbortAndExit(ERRORCODE_MPI_ERROR, "MPI_Send error");

	freeMatrix(&lenMatrix);
	freeMatrix(&result);
	MPI_Finalize();
	return 0;
}

