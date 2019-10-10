/** CITS3402 - High Performace Computing
 *  Project 2 - Shortest Paths
 *  Authors: Tomas Mijat 21721589 & Ethan Chin 22248878
 */

#include "paths.h"

    // Defining runCommand as global definition
#define runCommand "./paths"
#define SIZE 14


int main(int argc, char * argv[]) {

    FILE *fp;
    int *matrix;
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
    matrix = readFile(fp, &nelements, &dim);
    fclose(fp);

    // Printings
    printf("Dimensions: %d\n", dim);
    printf("Nelements: %d\n", nelements);

    for (int i = 0; i < nelements; i++)
    {
        if (i % dim == 0) printf("\n");
        printf("%d\t", matrix[i]);
    }
    printf("\n");

    free(matrix);
    return 0;
}

// Reads the file and allocates it to memory
int* readFile(FILE *fp, int *nelements, int *dim) {
    int num;
    int *matrix;

    // Reads the dimensions
    fread(&num, sizeof(int), 1, fp);
    *dim = num;

    // Allocate memory depending on dimensions
    matrix = calloc(sizeof(int), *dim * *dim);

    // Reads the elements
    while (fread(&num, sizeof(int), 1, fp) != 0) {
        matrix[(*nelements)++] = num;
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