COMPILATOR=mpicc
PROCESS_NUMBER=8
PROGRAM_NAME=merge

$(PROGRAM_NAME): mpi_$(PROGRAM_NAME).c ../lib/vector.c
	$(COMPILATOR) -o $@ mpi_$(PROGRAM_NAME).c ../lib/vector.c

clean:
	rm -f $(PROGRAM_NAME)
	rm -f mpi_$(PROGRAM_NAME).o

run: $(PROGRAM_NAME)
	mpirun -np $(PROCESS_NUMBER) $(PROGRAM_NAME) A.s rA.s

