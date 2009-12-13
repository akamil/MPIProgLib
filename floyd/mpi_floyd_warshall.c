#include <stdio.h>
#include <string.h>
#include <mpi.h>
#include <float.h>
#include "../lib/matrix.h"
#include "../lib/vector.h"

#define INFINITY  999999
#define PARSER_BUFFER 20
#define INFINITY_IN_FILE "-"

int ID = -1;
int PROCESS_NUMBER = -1;

#define TABLE_OF_CONSTANTS_SIZE 3
index_t tableOfConstants[TABLE_OF_CONSTANTS_SIZE] = { 0, 0, 0 };
#define NODES_NUMBER       tableOfConstants[0] /* graphs nodes numer */
#define X_MESH_SIZE        tableOfConstants[1] /* process number per all rows */
#define Y_MESH_SIZE        tableOfConstants[2] /* process number per all columns */

#define GRAPH_READFILE_NOERRORS       0
#define GRAPH_READFILE_CANT_OPEN_FILE 1
#define GRAPH_READFILE_CANT_MALLOC    2
#define GRAPH_READFILE_PARSE_ERROR    3
#define GRAPH_READFILE_MPI_ERROR      4

index_t calculateSubMatrixLength(index_t id, index_t mesh_size) {
	if (id < (mesh_size -1))
		return NODES_NUMBER / mesh_size;
	else
		return NODES_NUMBER - (id * (NODES_NUMBER / mesh_size));
}

inline index_t calculateSubMatrixRowLength(index_t id) {
	return calculateSubMatrixLength(id / X_MESH_SIZE, Y_MESH_SIZE);
}

inline index_t calculateSubMatrixColumnLength(index_t id) {
	return calculateSubMatrixLength(id % X_MESH_SIZE, X_MESH_SIZE);
}

index_t calculateSubMatrixBegin(index_t id, index_t mesh_size) {
	if (id == 0) return 0;
	if (id < (mesh_size -1))
		return (NODES_NUMBER / mesh_size) * (id -1);
	else
		return id * (NODES_NUMBER / mesh_size);
}

inline index_t calculateSubMatrixRowBegin(index_t id) {
	return calculateSubMatrixBegin(id / X_MESH_SIZE, Y_MESH_SIZE);
}

inline index_t calculateSubMatrixColumnBegin(index_t id) {
	return calculateSubMatrixBegin(id % X_MESH_SIZE, X_MESH_SIZE);
}

int BcastRowAndColumn(pvector_t row, pvector_t column, pmatrix_t subMatrix, index_t k) {
	static index_t indexOfProcessX = 0;
	static index_t indexOfProcessY = 0;
	index_t x, i, tmp;
	/* broadcast column k */
	if (k >= calculateSubMatrixColumnBegin(indexOfProcessX) + calculateSubMatrixColumnLength(indexOfProcessX))
		indexOfProcessX++;
	for(x = 0; x < Y_MESH_SIZE; x++) {
		tmp = indexOfProcessX + x * X_MESH_SIZE; /* number of next subprocess in column */
		if (tmp == ID)
			for(i = 0; i < subMatrix->rows; i++)
				setInVector(column, i + calculateSubMatrixRowBegin(ID), getFromMatrix(subMatrix, i, k - calculateSubMatrixColumnBegin(ID)));
		if (MPI_Bcast(&column->vector[calculateSubMatrixRowBegin(tmp)], calculateSubMatrixRowLength(tmp), MPI_UNIT, tmp, MPI_COMM_WORLD) != MPI_SUCCESS)
			return GRAPH_READFILE_MPI_ERROR;
	}
	/* broadcast row k */
	if (k >= calculateSubMatrixColumnBegin(indexOfProcessY) + calculateSubMatrixColumnLength(indexOfProcessY))
		indexOfProcessY++;
	for(x = 0; x < X_MESH_SIZE; x++) {
		tmp = indexOfProcessY * X_MESH_SIZE  + x; /* number of next subprocess in row */
		if (tmp == ID) {
			for(i = 0; i < subMatrix->columns; i++)
				setInVector(row, i + calculateSubMatrixColumnBegin(ID), getFromMatrix(subMatrix, k - calculateSubMatrixRowBegin(ID), i));
		}
		if (MPI_Bcast(&row->vector[calculateSubMatrixColumnBegin(tmp)], calculateSubMatrixColumnLength(tmp), MPI_UNIT, tmp, MPI_COMM_WORLD) != MPI_SUCCESS)
			return GRAPH_READFILE_MPI_ERROR;
	}
	return 0;
}

int Floyd_calculation(pmatrix_t lengthMatrix, index_t N) {
	index_t i, j, k, beginX, beginY;
	unit_t *tmp, element1, element2;
	pvector_t row, column;

	beginX = calculateSubMatrixColumnBegin(ID);
	beginY = calculateSubMatrixRowBegin(ID);

	initVector(&row, NODES_NUMBER);
	initVector(&column, NODES_NUMBER);
	if (!column || !row)
		return GRAPH_READFILE_CANT_MALLOC;
	for(k = 0; k < N; k++) {
		BcastRowAndColumn(row, column, lengthMatrix, k);
		for(i = 0; i < lengthMatrix->rows; i++)
			for(j = 0; j < lengthMatrix->columns; j++) {
				element1 = getFromVector(column, beginY + i);
				element2 = getFromVector(row, beginX + j);
				tmp = getMatrixElement(lengthMatrix, i, j);
				if ((element1 < INFINITY) && (element2 < INFINITY)) {
					element1 += element2;
					if (element1 < *tmp) *tmp = element1;
				}
			}
	}
	freeVector(&row);
	freeVector(&column);
	return 0;
}

void moveSubMatrixToMatrix(pmatrix_t source, pmatrix_t destiny) {
	index_t x, y;
	for(x = 0; x < destiny->columns; x++)
		for(y = 0; y < destiny->rows; y++)
			setInMatrix(destiny, y, x, getFromMatrix(source, y, x));
}

int readAndSendDataFromFile(FILE *file, index_t length_y, pmatrix_t subMatrix0, index_t IDInRow) {
	char          tmpRead[PARSER_BUFFER];
	index_t       x, y;
	pmatrix_t     tmpMatrix;
	MPI_Datatype  columntype;
	MPI_Status    stat;

	if (initMatrix(&tmpMatrix, length_y, NODES_NUMBER) == NULL)
		return GRAPH_READFILE_CANT_MALLOC;

	for(y = 0; y < length_y; y++)
		for(x = 0; x < NODES_NUMBER; x++) {
			if (fscanf(file, "%s", tmpRead) != 1)
				return GRAPH_READFILE_PARSE_ERROR;
			if (strcmp(tmpRead, INFINITY_IN_FILE) != 0)
				setInMatrix(tmpMatrix, y, x, atof(tmpRead));
			else
				setInMatrix(tmpMatrix, y, x, INFINITY);
		}
	length_y = 0;
	for(x = 0; x < X_MESH_SIZE; x++, length_y += calculateSubMatrixColumnLength(IDInRow), IDInRow++) {
		if (IDInRow == 0)
			moveSubMatrixToMatrix(tmpMatrix, subMatrix0);
		else {
			if (MPI_Type_vector(tmpMatrix->rows, calculateSubMatrixColumnLength(IDInRow), tmpMatrix->columns, MPI_UNIT, &columntype) != MPI_SUCCESS)
				return GRAPH_READFILE_MPI_ERROR;
			if (MPI_Type_commit(&columntype) != MPI_SUCCESS)
				return GRAPH_READFILE_MPI_ERROR;
			if (MPI_Send(&(tmpMatrix->tab[length_y]), 1, columntype, IDInRow, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
				return GRAPH_READFILE_MPI_ERROR;
		}
	}
	freeMatrix(&tmpMatrix);
	return GRAPH_READFILE_NOERRORS;
}

void printGraphLengthMatrix(pmatrix_t subMatrix, FILE *file) {
	index_t x, y;
	unit_t tmp;
	for(x = 0; x < subMatrix->rows; x++) {
		for(y = 0; y < subMatrix->columns; y++) {
			tmp = getFromMatrix(subMatrix, x, y);
			if (tmp == INFINITY)
				fprintf(file, INFINITY_IN_FILE" ");
			else
				fprintf(file, UNIT_SPECIFIER" ", tmp);
		}
		fprintf(file, "\n");
	}
}

int recivedAndWriteDataToFile(FILE *file, pmatrix_t subMatrix) {
	pmatrix_t tmpMatrix;
	index_t i, x;
	MPI_Datatype column;
	MPI_Status status;

	if (ID == 0) {
		fprintf(file, "%d\n", NODES_NUMBER);
		for(i = 0; i < PROCESS_NUMBER; i += X_MESH_SIZE) {
			initMatrix(&tmpMatrix, calculateSubMatrixRowLength(i), NODES_NUMBER);
			for(x = 0; x < X_MESH_SIZE; x++) {
				if (i + x == 0) {
				   	moveSubMatrixToMatrix(subMatrix, tmpMatrix);
				} else {
					if (MPI_Type_vector(calculateSubMatrixRowLength(i), calculateSubMatrixColumnLength(i+x), NODES_NUMBER, MPI_UNIT, &column) != MPI_SUCCESS)
						return GRAPH_READFILE_MPI_ERROR;
					if (MPI_Type_commit(&column) != MPI_SUCCESS)
						return GRAPH_READFILE_MPI_ERROR;
					if (MPI_Recv(&(tmpMatrix->tab[calculateSubMatrixColumnBegin(x)]), 1, column, i + x, 0, MPI_COMM_WORLD, &status) != MPI_SUCCESS)
						return GRAPH_READFILE_MPI_ERROR;
				}
			}
			printGraphLengthMatrix(tmpMatrix, file);
			freeMatrix(&tmpMatrix);
		}
	} else {
		if (MPI_Send(subMatrix->tab, getMatrixElementNumber(subMatrix), MPI_UNIT, 0, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
			return GRAPH_READFILE_MPI_ERROR;
	}
	return GRAPH_READFILE_NOERRORS;
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

#define ARGV_GRAPH_FILE_NAME     argv[1]
#define ARGV_MESH_X              argv[2]
#define ARGV_MESH_Y              argv[3]
#define ARGV_RESULT_FILE_NAME    argv[4]

int main(int argc, char** argv) {
	pmatrix_t lengthSubMatrix;
	FILE *file;
	int res;
  	index_t	i;
    MPI_Status status;

	if (argc < 4) {
		fprintf(stderr, "MPI Floyd-Warshall search for shortest paths.\nUsage:\n\t%s <graph file (to read)> <rows of processor mesh> <columns of processor mesh> <result graph file (to write)>\n", argv[0]);
		MPI_Finalize(); return 1;
	}

	if (mpiInit(&argc, &argv) != MPI_SUCCESS) {
		fprintf(stderr, "Can't start MPI\n");
		return 1;
	}

	if (ID == 0) {
		X_MESH_SIZE = atoi(ARGV_MESH_X);
		Y_MESH_SIZE = atoi(ARGV_MESH_Y);
		printf("Reading graph file: '%s'\n", ARGV_GRAPH_FILE_NAME);
		if ((file = fopen(ARGV_GRAPH_FILE_NAME, "r")) != NULL) {
			if (fscanf(file, "%d", &NODES_NUMBER) != 1)
				return AbortAndExit(GRAPH_READFILE_PARSE_ERROR, "Syntax error in file");
		} else
			return AbortAndExit(GRAPH_READFILE_CANT_OPEN_FILE, "Can't open graph file");
	}

	if ((res = MPI_Bcast(tableOfConstants, TABLE_OF_CONSTANTS_SIZE, MPI_INT, 0, MPI_COMM_WORLD)) != MPI_SUCCESS)
		return AbortAndExit(res, "MPI_BCast error");

	if (initMatrix(&lengthSubMatrix, calculateSubMatrixRowLength(ID), calculateSubMatrixColumnLength(ID)) == NULL)
		return AbortAndExit(GRAPH_READFILE_CANT_MALLOC, "Cannot allocate memory");

	if (ID == 0) {
		i = 0;
		while(i < PROCESS_NUMBER) {
			if ((res = readAndSendDataFromFile(file, calculateSubMatrixRowLength(i), lengthSubMatrix, i)) != GRAPH_READFILE_NOERRORS)
				return AbortAndExit(res, "File read error");
			i += X_MESH_SIZE;
		}
		fclose(file);
	} else {
		if ((res = MPI_Recv(lengthSubMatrix->tab, lengthSubMatrix->rows * lengthSubMatrix->columns, MPI_UNIT, 0, 0, MPI_COMM_WORLD, &status)) != MPI_SUCCESS)
			return AbortAndExit(res, "MPI_Recv error");
	}

	if ((res = Floyd_calculation(lengthSubMatrix, NODES_NUMBER)) != GRAPH_READFILE_NOERRORS)
		return AbortAndExit(res, "Error in calculation");

	if (ID == 0) {
		if ((file = fopen(ARGV_RESULT_FILE_NAME, "w")) != NULL) {
			if ((res = recivedAndWriteDataToFile(file, lengthSubMatrix)) != MPI_SUCCESS)
				AbortAndExit(res, "File write error");
			fclose(file);
		} else {
			fprintf(stderr, "Can't open result file. The results will be printed on standard output\n");
			if ((res = recivedAndWriteDataToFile(stdout, lengthSubMatrix)) != MPI_SUCCESS)
				AbortAndExit(res, "File write error");
		}
	} else if ((res = recivedAndWriteDataToFile(NULL, lengthSubMatrix)) != MPI_SUCCESS)
		AbortAndExit(res, "File write error");

	if (ID == 0) printf("Result file '%s' has been saved.\n", ARGV_RESULT_FILE_NAME);
	MPI_Finalize();
	return 0;
}

