/** CITS3402 - High Performace Computing
 *  Project 2 - Shortest Paths
 *  Authors: Tomas Mijat 21721589 & Ethan Chin 22248878
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"

// Function Declarations
int** readFile(FILE *fp, int *dim);
char* getFileName(int argCount, char *argInput[]);
void fileCheck(FILE *fp);