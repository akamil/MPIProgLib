#include "vector.h"

#define VECTOR_READFILE 
#define VECTOR_ERRORCODE_PARSER_ERROR 1

pvector_t initVector(pvector_t *vector, index_t length) {
	*vector = malloc(sizeof(vector_t));
	if (vector == NULL) return NULL;
	(*vector)->length = length;
	(*vector)->vector = calloc(length, sizeof(unit_t));
	if ((*vector)->vector == NULL) return NULL;
	return *vector;
}

inline pvector_t initCopyOfVector(pvector_t *matrix, const pvector_t original) {
	return initVector(matrix, original->length);
}

void freeVector(pvector_t *vector) {
	free((*vector)->vector);
	free(*vector);
}

inline unit_t *getVectorElement(pvector_t vector, index_t index) {
	return &(vector->vector[index]);
}

inline unit_t getFromVector(pvector_t vector, index_t index) {
	return *getVectorElement(vector, index);
}

inline void setInVector(pvector_t vector, index_t index, unit_t value) {
	*getVectorElement(vector, index) = value;
}

void printVectorOnStream(pvector_t vector, FILE *stream, const char *unit_t_format) {
	index_t i;
	fprintf(stream, "%d\n", vector->length);
	for(i = 0; i < vector->length; i++) 
		fprintf(stream, unit_t_format, getFromVector(vector, i));
	fprintf(stream, "\n");
}

inline void printVector(pvector_t vector) {
	printVectorOnStream(vector, stdout, UNIT_SPECIFIER " ");
}

int readVector(FILE *file, pvector_t *vector) {
	index_t number, x;
	if (fscanf(file, "%d\n", &number) == 1) {
		initVector(vector, number);
		for(x = 0; x < number; x++) {
			if (fscanf(file, UNIT_SPECIFIER "\n", getVectorElement(*vector, x)) != 1)
				return ERRORCODE_PARSER_ERROR;
		}
	} else
		return ERRORCODE_PARSER_ERROR;
	return ERRORCODE_NOERRORS;
}

