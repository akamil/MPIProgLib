COMPILATOR=mpicc
PROCES_NUMBER=4
PROGRAM_NAME=dijkstra

$(PROGRAM_NAME): mpi_$(PROGRAM_NAME).c ../lib/matrix.c ../lib/vector.c
	$(COMPILATOR) -o $@ mpi_$(PROGRAM_NAME).c ../lib/matrix.c ../lib/vector.c

clean:
	rm -f $(PROGRAM_NAME)
	rm -f mpi_$(PROGRAM_NAME).o

run: $(PROGRAM_NAME)
	mpirun -np $(PROCES_NUMBER) $(PROGRAM_NAME) D.g rD.g

