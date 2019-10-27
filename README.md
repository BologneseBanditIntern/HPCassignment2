# HPC Assignment2
### Ethan Chin - 22248878     Tomas Mijat - 21721589
-------------------------------------------------
#   How to run the Program

## Compiling the program

Within the folder there are a number of files. These files are;
- 'paths.c' : This is the source code for the program.
- 'paths.h' : This is the header file containing precompiler informationfor the program.
- 'makefile' : This is a tool used to easily compile the program. 
- 'examples' : This folder contains the exmpales provided by Nicolas Pritchard. 
 To compile the program use the command line instruction 'make' this will initiate the
compilation of the program, by running the following commands.
> mpicc -c paths.c
>
> mpicc -o paths paths.o

## Running the Program

To run the program use the command mpirun [-n {numberOfProcesses}] ./paths [-f] {fileName}
Where;  - objects within []'s are optional and may be included if wished.
        - nummberOfProcesses must equally divide the number of total datapoints. i.e. 5 nodes (25 data points) can be 1 or 5, 16 can be 1,2,4,8,16,32,etc
The program will print to the terminal timing data for each phase of execution (readingFile,ShortPathAlgorithm & writeFile) in addition to the total execution time, for each running process. This displace can be disabled by commenting out/deleting line 235 of the paths.c file and recompiling the program.
The program will also output a file containing the all pairs shortest path matrix (with the first digit representing the number of nodes of the graph), with the same name of the input file, with the exception that the extension is changed to ".out". This will be placed in the same file location as the input file. 


# Assignment Breif

The goal of this project is to design an implement parallel algorithms to solve
the all-pairs-shortest path problem for a number of large graphs using either
Dijkstra’s or Floyd-Warshall algorithms. You will use C with the MessagePassing Interface (MPI) to write parallelized code. You will also need to observe
and comment on how speedup is affected by changing the size of the problem and
the number of nodes you use. A reasonable range for the size of your problems
is 1024 to 4096 vertices. A reasonable range for the number of processors is one
to 16.

# Your submission will include:

# 1.1 Code
A copy of all source files, job scripts and build-files used. Source files must include your names and student numbers at the top. Your program must compile,
execute correctly and be well documented. Make sure to include comments to
explain any MPI function the first time it is used.

# 1.2 Report
Alongside your code you will also need to write a report addressing the following
points:
• How your approach partitions data and processing
• A table showing speedup for graphs of 2048 vertices vs 1, 2, 4, 8, 16 processors
• A table showing speedup for 4 processors vs graphs of 256, 512, 1024, 2048, 4096
vertices.
When discussing the performance of your approach, try to relate performance to
various factors in your computation (number of compute nodes, load balancing,
communication overhead etc.)

For further requirements please see:
http://teaching.csse.uwa.edu.au/units/CITS3402/labs/3402Assignment2.pdf
