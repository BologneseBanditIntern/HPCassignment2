/** CITS3402 - High Performace Computing
 *  Project 2 - Shortest Paths
 *  Authors: Tomas Mijat 21721589 & Ethan Chin 22248878
 */

#include "paths.h"

    // Defining runCommand as global definition
#define runCommand "./paths"
#define SIZE 14

#define MAX_ROWS 1000
#define MAX_COLS 1000


int main(int argc, char * argv[]) {

    FILE *fp;
    int **matrix;
    int *dist;
    int dim;
    int nelements = 0;

    //Conditional checking for invalid command line argument. Could potentially check for greater then 3 or 4 based on the -f flag, however this may be overkill.
    if (argc < 2 )
    {
        printf("An invalid number of arguments has be inputted. \nPlease check your command line arguments.\n");
        exit(0);
    }
    char * file = getFileName( argc, argv);
    printf("The file name is:\t%s\n%d\n", file,argc);

    // Open file
    fp = fopen(file, "rb");
    fileCheck(fp);

    // Reads the file data and sets the matrix array
    matrix = readFile(fp, &dim);
    fclose(fp); 

    // Initialize all distance balues as max (dim)
    nelements = dim*dim;
    dist = calloc(sizeof(int), nelements);

    // Prints the matrix
    printf("%d\n",dim);
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            printf("%d ",matrix[i][j]);
        }
        printf("\n");
    }

    free(matrix);
    return 0;
}

// Reads the file and allocates it to memory
int** readFile(FILE *fp, int *dim) {
    int num;
    int nelements = 0;
    int *ptr;
    int **matrix;

    // Reads the dimensions
    fread(&num, sizeof(int), 1, fp);
    *dim = num;

    // Allocate memory depending on dimensions
    matrix = (int **) malloc(sizeof(int *) * *dim + sizeof(int) * *dim * *dim);
    ptr = (int *)(matrix + *dim);

    // Adds row pointers
    for(int i = 0; i < *dim+1; i++) {
        matrix[i] = (ptr + *dim * i); 
    }

    // Reads the elements
    while (fread(&num, sizeof(int), 1, fp) != 0) {
        matrix[nelements / *dim][nelements % *dim] = num;
        nelements++;
    }
    
    return matrix;
}

//This function takes the command line arguments i.e. argc and argv and returns a char pointer, which is the filename/ file path.
char* getFileName(int argCount, char *argInput[])
{   
    int i;
    char * fileName;
    for(i = 0; i < argCount;i++){       //For loop finding the file name argument utilizing both -f flags and without.
        if(strcmp(argInput[i],runCommand)==0) 
        {
            if(strcmp(argInput[i+1],"-f")==0)
            {
                fileName = argInput[i+2];
            }
            else
            {
                fileName = argInput[i+1];   
            }
        }
    }
    return fileName;
}

// Checks that the file can be opened
void fileCheck(FILE *fp) {
    if (fp == NULL) {
        perror("Error: This file cannot be accessed.");
        exit(EXIT_FAILURE);
    }
    return;
}



/** MPI Broadcast
         *  One process sends the same information to every other process,
         *  OpenMPI chooses the most optimal algorithm depending on the conditions
         *  MPI_Bcast(void *buffer, int count, MPI_Datatype, int root, MPI_Comm comm)
         *  @param buffer starting address of buffer
         *  @param count number of entries in buffer
         *  @param MPI_Datatype data type of the buffer
         *  @param root rank of broadcast root
         *  @param comm communicator
        */
        //mpierror = MPI_Bcast(b, dim, MPI_DOUBLE, root, MPI_COMM_WORLD);