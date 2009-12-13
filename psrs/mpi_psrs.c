#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "../lib/vector.h"

#define ROOT_ID 0
int PROCESS_NUMBER = -1;
int ID = -1;

#define TABLE_OF_CONSTANTS_SIZE 3
index_t tableOfConstants[TABLE_OF_CONSTANTS_SIZE] = { 0, 0, 0 };
#define ELEMENTS_NUMBER      tableOfConstants[0] /* size of data set */
#define ELEMENTS_PER_PROCESS tableOfConstants[1] /* size of this processor's data subset */
#define SQR_PROCESS_NUMBER   tableOfConstants[2] /* square of process number */

int initMPI(int *argc, char ***argv) {
	int res;
	if ((res = MPI_Init(argc, argv)) != MPI_SUCCESS) return res;
	if ((res = MPI_Comm_rank(MPI_COMM_WORLD, &ID)) != MPI_SUCCESS) return res;
	if ((res = MPI_Comm_size(MPI_COMM_WORLD, &PROCESS_NUMBER)) != MPI_SUCCESS) return res;
	SQR_PROCESS_NUMBER = PROCESS_NUMBER * PROCESS_NUMBER;
	return MPI_SUCCESS;
}

index_t listLength(index_t id) {
	if (id < (PROCESS_NUMBER -1))
		return (ELEMENTS_NUMBER / PROCESS_NUMBER);
	else
		return (ELEMENTS_NUMBER - (ELEMENTS_NUMBER/PROCESS_NUMBER) * (PROCESS_NUMBER -1));
}

void quicksort(unit_t *T, int low, int high) {
	int i, j;
	unit_t x, tmp;

	x = T[(low + high)/2];
	i = low;
	j = high;
	do {
		while (T[i] < x) ++i;
		while (T[j] > x) --j;
		if (i<=j) {
			tmp = T[i];
			T[i] = T[j];
			T[j] = tmp;
			++i; --j;
		}
	} while(i < j);
	if (low < j) quicksort(T, low, j);
	if (high > i) quicksort(T, i, high);
}

inline void quicksortVector(pvector_t vector) {
	quicksort(vector->vector, 0, vector->length - 1);
}

void copyVector(pvector_t source, pvector_t destiny, index_t length) {
	index_t i;
	for(i = 0; i < length; i++)
		setInVector(destiny, i, getFromVector(source, i));
}

void copySampleToVector(pvector_t source, pvector_t destiny, index_t by, index_t number) {
	int step, i;
	for(i = 0, step = by; i < number; i++, step += by)
		setInVector(destiny, i, getFromVector(source, step));
}

int dataExchange(unit_t tmin, unit_t tmax, pvector_t *v, pvector_t tmp) {
	index_t i, x, pos = 0;
	unit_t t;
	pvector_t data;

	if (initVector(&data, ELEMENTS_PER_PROCESS << 1) == NULL)
		return ERRORCODE_CANT_MALLOC;

	for(i = 0; i < PROCESS_NUMBER; i++) {
		if (ID == i) copyVector(*v, tmp, (*v)->length);

		if (MPI_Bcast(tmp->vector, listLength(i), MPI_UNIT, i, MPI_COMM_WORLD) != MPI_SUCCESS)
			return ERRORCODE_MPI_ERROR;

		for(x = 0; x < listLength(i); x++) {
			t = getFromVector(tmp, x);
			if (t <= tmax) {
				if (t > tmin)
					setInVector(data, pos++, t);
			} else
				x+= ELEMENTS_NUMBER;
		}
	}

	freeVector(v);
	if (initVector(v, pos) == NULL)
		return ERRORCODE_CANT_MALLOC;
	copyVector(data, *v, pos);
	freeVector(&data);
	return ERRORCODE_NOERRORS;
}

pvector_t openVectorFile(char *fileName) {
	pvector_t v;
	FILE *file;

	printf("Reading set file '%s'\n", fileName);
	if ((file = fopen(fileName, "r")) != NULL) {
		if (readVector(file, &v) != ERRORCODE_NOERRORS)
			v = NULL;
		fclose(file);
	}
	return v;
}

void writeVectorToFile(pvector_t vector, char *fileName) {
	FILE *file;

	if ((file = fopen(fileName, "w")) != NULL) {
		printVectorOnStream(vector, file, UNIT_SPECIFIER "\n");
		printf("Result file '%s' has been saved.\n", fileName);
	} else {
		fprintf(stderr, "Can't open file '%s' to write, outcome data will be written on standard output\n", fileName);
		printVector(vector);
	}
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
	pvector_t v, tmp = NULL, samples = NULL;
	index_t i, length, step;
	unit_t min, max;
	MPI_Status status;
	MPI_Datatype sampleDatatype;

	if (initMPI(&argc, &argv) != MPI_SUCCESS)
		return AbortAndExit(ERRORCODE_MPI_ERROR, "Cannot initialize MPI.");

	if (argc < 3) {
		fprintf(stderr, "MPI Parallel Sorting by Regular Sampling implementation.\nUsage:\n\t%s <data set (to read)> <result  file (to write)>\n", argv[0]);
		MPI_Finalize(); return 1;
	}

	if (ID == ROOT_ID) {
		tmp = openVectorFile(ARGV_FILE_NAME);
		printf("Data set size: %d, process number: %d\n", tmp->length, PROCESS_NUMBER);
		if ((tmp->length/PROCESS_NUMBER) <= PROCESS_NUMBER)
			AbortAndExit(ERRORCODE_SIZE_DONT_MATCH, "Processor number is too big or size of data set is too small for correct calculation.\n");
		ELEMENTS_NUMBER = tmp->length;
	}

	if (MPI_Bcast(tableOfConstants, TABLE_OF_CONSTANTS_SIZE, MPI_INT, ROOT_ID, MPI_COMM_WORLD) != MPI_SUCCESS)
		return AbortAndExit(ERRORCODE_MPI_ERROR, "MPI_Bcast error.");

	ELEMENTS_PER_PROCESS = listLength(ID);
	initVector(&v, ELEMENTS_PER_PROCESS);

	if (ID == ROOT_ID) { /* Bcast data set */
		copyVector(tmp, v, v->length);
		for(i = 1, step = ELEMENTS_PER_PROCESS; i < PROCESS_NUMBER; i++) {
			if (MPI_Send(&(tmp->vector[step]), listLength(i), MPI_UNIT, i, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
				return AbortAndExit(ERRORCODE_MPI_ERROR, "MPI_Send error.");
			step += listLength(i);
		}
	} else if (MPI_Recv(v->vector, ELEMENTS_PER_PROCESS, MPI_UNIT, ROOT_ID, 0, MPI_COMM_WORLD, &status) != MPI_SUCCESS)
		return AbortAndExit(ERRORCODE_MPI_ERROR, "MPI_Recv error.");

	quicksortVector(v);

	if (initVector(&samples, PROCESS_NUMBER -1) == NULL)
		return AbortAndExit(ERRORCODE_CANT_MALLOC, "Cannot allocate memory for samples vector.");

	MPI_Type_vector(PROCESS_NUMBER, 1, ELEMENTS_NUMBER / SQR_PROCESS_NUMBER, MPI_UNIT, &sampleDatatype);
	MPI_Type_commit(&sampleDatatype);

	if (ID != ROOT_ID) { /* Sending samples to root proces */
 		if (MPI_Send(v->vector, 1, sampleDatatype, ROOT_ID, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
			return AbortAndExit(ERRORCODE_MPI_ERROR, "MPI_Send error.");
		if (initVector(&tmp, listLength(PROCESS_NUMBER -1)) == NULL)
			return AbortAndExit(ERRORCODE_CANT_MALLOC, "Cannot allocate memory for temporary vector.");
	} else { /* Reciving samples */
		copySampleToVector(v, tmp, (v->length)/PROCESS_NUMBER, PROCESS_NUMBER);
		for(step = PROCESS_NUMBER, i = 1; i < PROCESS_NUMBER; i++, step += PROCESS_NUMBER)
			if (MPI_Recv(&(tmp->vector[step]), PROCESS_NUMBER, MPI_UNIT, i, 0, MPI_COMM_WORLD, &status) != MPI_SUCCESS)
				return AbortAndExit(ERRORCODE_MPI_ERROR, "MPI_Recv error.");
		quicksort(tmp->vector, 0, SQR_PROCESS_NUMBER);
		copySampleToVector(tmp, samples, SQR_PROCESS_NUMBER / (PROCESS_NUMBER - 1), PROCESS_NUMBER - 1);
	}

	/* Broadcast selected samples to processors */
	if (MPI_Bcast(samples->vector, PROCESS_NUMBER-1, MPI_UNIT, ROOT_ID, MPI_COMM_WORLD) != MPI_SUCCESS)
		return AbortAndExit(ERRORCODE_MPI_ERROR, "MPI_Bcast error.");

	if ((i = dataExchange((ID == 0) ? UNITT_MIN : getFromVector(samples, ID -1), (ID == (PROCESS_NUMBER - 1)) ? UNITT_MAX : getFromVector(samples, ID), &v, tmp)) != ERRORCODE_NOERRORS)
		return AbortAndExit(i, "Error in while of data exchange.");

	/* Sorting new data */
	quicksortVector(v);

	if (ID != ROOT_ID) { /* Sending sorted data */
		if (MPI_Send(&(v->length), 1, MPI_INT, ROOT_ID, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
			return AbortAndExit(ERRORCODE_MPI_ERROR, "MPI_Send (sending size of data) error.");
		if (MPI_Send(v->vector, v->length, MPI_UNIT, ROOT_ID, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
			return AbortAndExit(ERRORCODE_MPI_ERROR, "MPI_Send error.");
	} else { /* Receiving sorted data */
		copyVector(v, tmp, v->length);
		for(step = v->length, i = 1; i < PROCESS_NUMBER; i++) {
			if (MPI_Recv(&length, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &status) != MPI_SUCCESS)
				return AbortAndExit(ERRORCODE_MPI_ERROR, "MPI_Recv (sending size of data) error.");
			if (MPI_Recv(&(tmp->vector[step]), length, MPI_UNIT, i, 0, MPI_COMM_WORLD, &status) != MPI_SUCCESS)
				return AbortAndExit(ERRORCODE_MPI_ERROR, "MPI_Recv error.");
			step += length;
		}
		writeVectorToFile(tmp, ARGV_RESULT_NAME);
		freeVector(&tmp);
	}
	freeVector(&v);
	MPI_Finalize();
	return 0;
}
