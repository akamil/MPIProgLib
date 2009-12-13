#include <stdio.h>
#include <mpi.h>
#include "../lib/vector.h"

int PROCESS_NUMBER = -1;
int ID = -1;

#define TABLE_OF_CONSTANTS_SIZE 3
index_t tableOfConstants[TABLE_OF_CONSTANTS_SIZE] = { 0, 0, 0 };
#define ELEMENTS_NUMBER      tableOfConstants[0] /* size of data set */
#define USE_PROCESS_NUMER    tableOfConstants[1] /* number of used processors */
#define ELEMENTS_PER_PROCESS tableOfConstants[2] /* size od subset for this processor */

void singleMerge(unit_t *a, index_t low, index_t high, index_t mid) {
	int i, j, k, c[50];
	i= low;
	j= mid + 1;
	k= low;

	while((i <= mid) && (j <= high)) {
		if (a[i] < a[j]) {
			c[k] = a[i];
			k++;
			i++;
		} else {
			c[k] = a[j];
			k++;
			j++;
		}
	}
	while(i <= mid) {
		c[k] = a[i];
		k++;
		i++;
	}
	while(j <= high) {
		c[k] = a[j];
		k++;
		j++;
	}
	for(i=low;i<k;i++)
		a[i] = c[i];
}

void mergesort(unit_t *tab, index_t low, index_t high) {
	int mid;
	if (low < high) {
		mid = (low+high)/2;
		mergesort(tab, low, mid);
		mergesort(tab, mid+1, high);
		singleMerge(tab, low, high, mid);
	}
	return;
}

pvector_t merge(pvector_t left, index_t lengthOfLeft, pvector_t right) {
	pvector_t tmp;
	index_t ileft = 0, iright = 0, itmp = 0;
	unit_t *p;

	if (initVector(&tmp, lengthOfLeft + right->length) == NULL)
		return NULL;

	while(ileft < lengthOfLeft && iright < right->length) {
		if (getFromVector(left, ileft) < getFromVector(right, iright))
			setInVector(tmp, itmp, getFromVector(left, ileft++));
		else
			setInVector(tmp, itmp, getFromVector(right, iright++));
		itmp++;
	}

	if (ileft == lengthOfLeft)
		p = getVectorElement(right, iright);
	else
		p = getVectorElement(left, ileft);

	for(; itmp < tmp->length; itmp++, p++)
		setInVector(tmp, itmp, *p);

	freeVector(&left);
	return tmp;
}

int openVectorFile(char *fileName, pvector_t *vector) {
	int res;
	FILE *file;

	printf("Reading set file '%s'\n", fileName);
	if ((file = fopen(fileName, "r")) != NULL) {
		res = readVector(file, vector);
		fclose(file);
		return res;
	} else
		return ERRORCODE_CANT_FILE_OPEN;
}

void writeVectorToFile(pvector_t vector, char *fileName) {
	FILE *file;

	if ((file = fopen(fileName, "w")) != NULL) {
		printVectorOnStream(vector, file, UNIT_SPECIFIER "\n");
		printf("Result file '%s' has been saved.\n", fileName);
	} else {
		fprintf(stderr, "Can't open file '%s' to write, matrix will be written on standard output\n", fileName);
		printVector(vector);
	}
}

int AbortAndExit(int result, char *reason) {
	fprintf(stderr, "%s\n", reason);
	MPI_Abort(MPI_COMM_WORLD, result);
	return result;
}

index_t listLength(index_t id) {
	if (id < (USE_PROCESS_NUMER -1))
		return (ELEMENTS_NUMBER / USE_PROCESS_NUMER);
	else
		return (ELEMENTS_NUMBER - (ELEMENTS_NUMBER/USE_PROCESS_NUMER) * (USE_PROCESS_NUMER -1));
}

index_t calculateListLength(index_t startID, index_t level) {
	index_t i, result;
	for(result = 0, i = startID; i < startID + level; i++)
		result += listLength(i);
	return result;
}

index_t sendTo(index_t level) {
	int i = ID;
	while(i % level != 0) i--;
	return i;
}

int initMPI(int *argc, char ***argv) {
	int res, s2;
	if ((res = MPI_Init(argc, argv)) != MPI_SUCCESS) return res;
	if ((res = MPI_Comm_rank(MPI_COMM_WORLD, &ID)) != MPI_SUCCESS) return res;
	if ((res = MPI_Comm_size(MPI_COMM_WORLD, &PROCESS_NUMBER)) != MPI_SUCCESS) return res;

	s2 = 1;
	while(s2 <= PROCESS_NUMBER) s2 *= 2;
	USE_PROCESS_NUMER = s2/2;
	return MPI_SUCCESS;
}

#define ARGV_PROGRAM_NAME argv[0]
#define ARGV_FILE_NAME    argv[1]
#define ARGV_RESULT_NAME  argv[2]

int main(int argc, char **argv) {
	pvector_t v, tmp;
	index_t i, b, level, prevLevel, len;
	MPI_Status status;

	if (initMPI(&argc, &argv) != MPI_SUCCESS)
		return AbortAndExit(ERRORCODE_MPI_ERROR, "Cannot initialize MPI.");

	if (argc < 3) {
		fprintf(stderr, "MPI Parallel Mergesort algorithm implementation.\nUsage:\n\t%s <data set (to read)> <result  file (to write)>\n", argv[0]);
		MPI_Finalize(); return 1;
	}

	if (ID == 0) {
		if (openVectorFile(ARGV_FILE_NAME, &v) != ERRORCODE_NOERRORS)
			AbortAndExit(ERRORCODE_CANT_READ_DATA, "Cannot read data file.");
		ELEMENTS_NUMBER = v->length;
	}

	if (MPI_Bcast(tableOfConstants, TABLE_OF_CONSTANTS_SIZE, MPI_INT, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
		return AbortAndExit(ERRORCODE_MPI_ERROR, "MPI_Bcast error.");

	ELEMENTS_PER_PROCESS = listLength(ID);

	if (ID == 0) { /* Bcast data */
		for(i = 1, b = ELEMENTS_PER_PROCESS; i < USE_PROCESS_NUMER; i++) {
			if (MPI_Send(&(v->vector[b]), listLength(i), MPI_UNIT, i, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
				return AbortAndExit(ERRORCODE_MPI_ERROR, "MPI_Send error.");
			b += listLength(i);
		}
	} else { /* Recive data from main proces */
		initVector(&v, ELEMENTS_PER_PROCESS);
		if (MPI_Recv(v->vector, ELEMENTS_PER_PROCESS, MPI_UNIT, 0, 0, MPI_COMM_WORLD, &status) != MPI_SUCCESS)
			return AbortAndExit(ERRORCODE_MPI_ERROR, "MPI_Recv error.");
	}

	/* Major calculation */
	mergesort(v->vector, 0, ELEMENTS_PER_PROCESS -1);

	/* Merge phase */
	for(prevLevel = 1, level = 2; level <= USE_PROCESS_NUMER; level <<= 1, prevLevel <<= 1) {
		if (ID % level == 0) {
			if (initVector(&tmp, calculateListLength(ID + prevLevel, prevLevel)) == NULL)
				return AbortAndExit(ERRORCODE_CANT_MALLOC, "Cannot allocate memory for temporary vector.");
			if (MPI_Recv(tmp->vector, tmp->length, MPI_UNIT, ID + prevLevel, 0, MPI_COMM_WORLD, &status) != MPI_SUCCESS)
				return AbortAndExit(ERRORCODE_MPI_ERROR, "MPI_Recv error.");
		} else {
			if (MPI_Send(v->vector, v->length, MPI_UNIT, sendTo(level), 0, MPI_COMM_WORLD) != MPI_SUCCESS)
				return AbortAndExit(ERRORCODE_MPI_ERROR, "MPI_Send error.");
			freeVector(&v); /* This process don't need vector v any more */
			break; /* Becase, this process is ending now */
		}
		v = merge(v, (level == 2) ? ELEMENTS_PER_PROCESS : v->length, tmp);
		if (v == NULL)
			return AbortAndExit(ERRORCODE_CANT_MALLOC, "Cannot allocate memory in merge function.");
	}

	if (ID == 0) {
		writeVectorToFile(v, ARGV_RESULT_NAME);
		freeVector(&v);
	}

	MPI_Finalize();
	return 0;
}
