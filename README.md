# HPC Assignment2

-------------------------------------------------

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
