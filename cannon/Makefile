COMPILATOR=mpicc
PROCES_NUMBER=4

cannon: mpi_cannon.c ../lib/matrix.c
	$(COMPILATOR) -o $@ mpi_cannon.c ../lib/matrix.c

clean:
	rm -f cannon
	rm -f mpi_cannon.o

run: cannon
	mpirun -np $(PROCES_NUMBER) cannon A.m B.m AB.m

