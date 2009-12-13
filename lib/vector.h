#ifndef __VECTOR_H__
#define __VECTOR_H__

#include <stdlib.h>
#include <stdio.h>
#include "unit.h"
#include "error.h"

typedef struct vector_t {
	unit_t   *vector;
	index_t   length;
} vector_t;

typedef vector_t *pvector_t;

pvector_t       initVector(pvector_t *vector, index_t length);
pvector_t       initCopyOfVector(pvector_t *matrix, const pvector_t original);
void            freeVector(pvector_t *vector);

inline unit_t  *getVectorElement(pvector_t vector, index_t index);
inline unit_t   getFromVector(pvector_t vector, index_t index);
inline void     setInVector(pvector_t vector, index_t index, unit_t value);

void            printVectorOnStream(pvector_t vector, FILE * stream, const char *unit_t_format);
inline void     printVector(pvector_t vector);

int             readVector(FILE *file, pvector_t *vector);

#endif

