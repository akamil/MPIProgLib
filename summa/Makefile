COMPILATOR=mpicc
PROCES_NUMBER=9
PROGRAM_NAME=summa

$(PROGRAM_NAME): mpi_summa.c ../lib/matrix.c
	$(COMPILATOR) -o $@ mpi_summa.c ../lib/matrix.c

clean:
	rm -f $(PROGRAM_NAME)
	rm -f mpi_suma.o

run: $(PROGRAM_NAME)
	mpirun -np $(PROCES_NUMBER) $(PROGRAM_NAME) A.m B.m AB.m

