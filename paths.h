/** CITS3402 - High Performace Computing
 *  Project 2 - Shortest Paths
 *  Authors: Tomas Mijat 21721589 & Ethan Chin 22248878
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "mpi.h"

// Function Declarations
int* readFile(FILE *fp, int *dim);
int readFileDims(FILE *fp);
char* getFileName(int argCount, char *argInput[]);
int* dijkstraP(int *matrix, int dim, int local_n, int myRank, int *root_matrix);
int* dijkstra(int *matrix, int dim);
int* initMatrix(int dim);
int* initMatrixP(int dim);
void memory_check(int *matrix, char *msg);
void fileCheck(FILE *fp);
void mpi_error_check(int mpierror);
