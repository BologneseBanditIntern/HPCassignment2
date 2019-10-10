# CITS3402 - High Performace Computing
# Project 2 - Shortest Paths
# Authors: Tomas Mijat 21721589 & Ethan Chin 22248878

OBJ     =   paths.o
PROJECT =   paths
mpi = mpicc

$(PROJECT) : $(OBJ)
	$(mpi) -o $(PROJECT) $(OBJ)
	@echo "Use 'mpirun ./paths <arg>' to run."

paths.o : paths.c
	$(mpi) -c paths.c

clean:
	@echo "Cleaning up..."
	rm *.o
	rm *.exe