/** CITS3402 - High Performace Computing
 *  Project 2 - Shortest Paths
 *  Authors: Tomas Mijat 21721589 & Ethan Chin 22248878
 */

#include "paths.h"

    // Defining runCommand as global definition
#define runCommand "./paths"


int main(int argc, char * argv[]) {
    printf("Hello!\n");
        //Conditional checking for invalid command line argument. Could potentially check for greater then 3 or 4 based on the -f flag, however this may be overkill.
    if (argc < 2 )
    {
        printf("An invalid number of arguments has be inputted. \nPlease check your command line arguments.\n");
        exit(0);
    }
    char * file = getFileName( argc, argv);
    printf("The file name is:\t%s\n%d\n", file,argc);
    return 0;
}



char* getFileName(int argCount, char *argInput[])
{   //This function takes the command line arguments i.e. argc and argv and returns a char pointer, which is the filename/ file path.
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