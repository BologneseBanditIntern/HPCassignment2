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
    MPI_File fileHandle;        //MPI File type
    int *root_matrix;
    int *root_dist;
    int *local_matrix;
    int *local_dist;
    int dim;
    int nelements = 0;
    int myRank, numProcs, mpierror;
    int root = 0;
    int local_n;
    int *dist;
    clock_t shortPath_start, shortPath_end, file_start, file_end;

    file_start = clock();

    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    //Conditional checking for invalid command line argument. Could potentially check for greater then 3 or 4 based on the -f flag, however this may be overkill.
    if (argc < 2 )
    {
        printf("An invalid number of arguments has be inputted. \nPlease check your command line arguments.\n");
        exit(0);
    }

    char * file = getFileName( argc, argv );
    if (myRank == root) printf("The file name is:\t%s\n", file);

    if(myRank == root)  //if root process, then find the number of dimensions
    {
        FILE *fileDims = fopen(file,"rb");
        if(fileDims == NULL)
        {
            printf("The inputted file does not exist or cannot be openned. Please check inputs and  try again.\n");
            dim = -1;
        }
        else
        {
        dim = readFileDims(fileDims);
        fclose(fileDims);
        }
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
    // Broadcast the numebr of dimensions
    mpierror = MPI_Bcast(&dim, 1, MPI_INT, root, MPI_COMM_WORLD); mpi_error_check(mpierror);
    MPI_Barrier(MPI_COMM_WORLD);
    if(dim == -1){      //If file can not be accessed, gracefully terminate.
        MPI_Finalize();
        return 0;
    }
    nelements = dim * dim;

    // Open file
    MPI_File_open(MPI_COMM_WORLD,file,MPI_MODE_RDONLY,MPI_INFO_NULL,&fileHandle);

    //Use collective call to read file
    int elementsPerProcess = nelements / numProcs;
    while((elementsPerProcess * numProcs) < nelements) //Dividing work evenly (ish) Idea is that if the remaining parts of data is too great, then increase the amount of data each process should gather by 1
    {
        elementsPerProcess++;
    }
    int elementsLastProcess = nelements - (elementsPerProcess * (numProcs-1));      //Assigning a smaller number to the last process

    //if(myRank == root)  printf("The number of elements read per process is:\t%d\nWith the last process reading:\t %d  elements\n",elementsPerProcess,elementsLastProcess);

    //Reading of the file
    int * partialMatrix;


    //Allocating memory to the for partial matrix then reading file to processes
    if(myRank == (numProcs - 1) && numProcs > 1)
    {
        MPI_Offset offset = 1 + (myRank * elementsPerProcess);
        partialMatrix = malloc(elementsLastProcess*sizeof(int));
        MPI_File_sync(fileHandle);
        mpierror = MPI_File_read_at(fileHandle,offset*sizeof(int),partialMatrix,elementsLastProcess,MPI_INT,&status);
        mpi_error_check(mpierror);
    }
    else
    {
        MPI_Offset offset = 1 + (myRank * elementsPerProcess);
        partialMatrix = malloc(elementsPerProcess * sizeof(int));
        MPI_File_sync(fileHandle);
        mpierror = MPI_File_read_at(fileHandle,offset*sizeof(int),partialMatrix,elementsPerProcess,MPI_INT,&status);
        mpi_error_check(mpierror);
    }

    /*
    if(mpierror == MPI_SUCCESS)
    {
        printf("\nThe read was successful!\n\n");
    }
    */

    MPI_File_close(&fileHandle);        //Closing the file after finishing reading.

    //Giving all processes the contents of the file
    root_matrix = (int *) malloc((dim*dim) * sizeof(int));     //Allocated pointer space for root matrix.

    if(myRank == (numProcs - 1) && numProcs > 1)
    {   //All gather statement for the final process (might have less elements)
        mpierror = MPI_Allgather(partialMatrix,elementsLastProcess,MPI_INT,root_matrix,elementsPerProcess,MPI_INT,MPI_COMM_WORLD);
        mpi_error_check(mpierror);
    }
    else
    {   //All gather statement for every other process
        mpierror = MPI_Allgather(partialMatrix,elementsPerProcess,MPI_INT,root_matrix,elementsPerProcess,MPI_INT,MPI_COMM_WORLD);
        mpi_error_check(mpierror);
         
    }

    file_end = clock();
    double timespentFile = (double) (file_end - file_start) / CLOCKS_PER_SEC;       //Measure time spent reading the input file

    if (myRank == root) {
        /*
        printf("File matrix:\n");
        for (int i = 0; i < nelements; i++)
        {
            printf("%d ", root_matrix[i]);
            if (i % dim == dim-1) printf("\n");
        }
        */
        
        root_dist = initMatrix(dim);
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
    // broadcasts the number of dimensions
    //mpierror = MPI_Bcast(&dim, 1, MPI_INT, root, MPI_COMM_WORLD); mpi_error_check(mpierror);
    //printf("DIM: %d %d\n", myRank, dim);

    // Each process will get a piece of the array
    local_n = (dim * dim) / numProcs; // p / n
    //if (myRank == root) printf("local:%d\n", local_n);

    // Start time of shortest path operations
    shortPath_start = clock();

    if (myRank == root) {
        //root_dist = dijkstra(root_matrix, dim);
    }

    local_dist = dijkstraP(dim, local_n, myRank, root_matrix);

    /** MPI Gather
     *  Collects all the results from the processes into the root process.
     *  Processors send their elements to the root array to be collected.
     *  The elements are ordered by the rank of the process
     *  local_n is the number of elements received per process
     *  local_dist is the elements stored on the process which collected to the root_dist
     */ 
    mpierror = MPI_Gather(local_dist, local_n, MPI_INT, root_dist, local_n, MPI_INT, root, MPI_COMM_WORLD);
    mpi_error_check(mpierror);

    /*
    if (myRank == root) {
        // Prints the distance matrix
        printDistance(root_dist, dim, nelements);
    }
    */
    

    // Prints the time of the execution
    shortPath_end = clock();
    double timespentShortPath = (double) (shortPath_end - shortPath_start) / CLOCKS_PER_SEC;
    printf("***Process %d****\nTime spent on file reading:\t%lf\nTime spent on shortest path op:\t%lf\n\n",myRank,timespentFile,timespentShortPath);
    
    if(myRank == root) free(root_dist);
    free(root_matrix);
    free(local_dist);
    

    
    MPI_Finalize();


    return 0;
}

// Initializes a 1d matrix size of dim x dim
int* initMatrix(int dim) {
    int *matrix;

    // Allocate memory depending on dimensions
    matrix = (int *) malloc(sizeof(int) * dim * dim);
    memory_check(matrix, "initializing matrix");
    return matrix;
}

// Initializes a 1d matrix size of dim
int* initMatrixP(int dim) {
    int *matrix;

    // Allocate memory depending on dimensions
    matrix = (int *) malloc(sizeof(int) * dim);
    memory_check(matrix, "initializing parallel matrix");
    return matrix;
}

int* dijkstraP(int dim, int local_n, int myRank, int *root_matrix) {

    // output array, holds the shortest distances
    int *dist = initMatrixP(local_n);

    // 1 if shortest distance has been found
    int *visited = initMatrixP(local_n);

    int pos;
    int rows = local_n / dim;

    // Initialize all distance values as max (dim) and visited
    for (int i = 0; i < local_n; i++)
    {
        dist[i] = dim;
        visited[i] = 0;
    }
    //printf("My Rank: %d:  %d, %d, %d, %d\n", myRank, matrix[0], matrix[1], matrix[2] ,matrix[3]);
    // iterates through all vertices
    for (int n = 0; n < rows; n++) {

        // distance from self is 0
        dist[(n * dim + n) + myRank * rows] = 0;

        // Finds the shortest paths for all verticies from n
        for (int count = 0; count < dim - 1; count++)
        {
            // Finds the min dist value, can be moved to a separate function
            //int u = minDistance(dist, visited);
            int min = dim;
            int u;

            for (int v = 0; v < dim; v++) {
                pos = n * dim + v;
                if (visited[pos] == 0 && dist[pos] <= min) {
                    min = dist[pos];
                    u = v;
                }
            }
            // end of min

            // set vertex as visited
            visited[n * dim + u] = 1;
            //printf("U: %d %d\n", u, myRank);

            for (int v = 0; v < dim; v++)
            {
                pos = n * dim + v;
                if (!visited[pos] && root_matrix[u * dim + v] && dist[n * dim + u] != dim &&
                    dist[n * dim + u] + root_matrix[u * dim + v] < dist[pos]) {
                    //printf("%d, %d, %d, %d, %d, %d, %d, %d\n", visited[pos], root_matrix[u * dim + v], dist[n*dim+u], dist[pos], v , u, pos, n);
                    dist[pos] = dist[n * dim + u] + root_matrix[u * dim + v];
                }
            }
        }
        // print
        /*
        printf("Vertiex\t Distance\n");
        for (int i = 0; i < dim; i++)
        {
            printf("%d ", dist[i]);
        }
        */
        
    }
    free(visited);

    return dist;
}


int* dijkstra(int *matrix, int dim) {

    // output array, holds the shortest distances
    int *dist = initMatrix(dim);

    // 1 if shortest distance has been found
    int *visited = initMatrix(dim);

    int pos;

    // Initialize all distance values as max (dim) and visited
    for (int i = 0; i < dim; i++)
    {
        for (int j = 0; j < dim; j++)
        {
            pos = i * dim + j;
            dist[pos] = dim;
            visited[pos] = 0;
        }
    }
    
    // iterates through all vertices
    for (int n = 0; n < dim; n++) {

        // distance from self is 0
        dist[n * dim + n] = 0;

        // Finds the shortest paths for all verticies from n
        for (int count = 0; count < dim - 1; count++)
        {
            // Finds the min dist value, can be moved to a separate function
            //int u = minDistance(dist, visited);
            int min = dim;
            int min_index;

            for (int v = 0; v < dim; v++) {
                pos = n * dim + v;
                if (visited[pos] == 0 && dist[pos] <= min) {
                    min = dist[pos];
                    min_index = v;
                }
            }
            // end of min

            // set vertex as visited
            int u = min_index;
            visited[n * dim + u] = 1;
            //printf("U: %d\n", u);

            for (int v = 0; v < dim; v++)
            {
                pos = n * dim + v;
                printf("%d, %d, %d, %d, %d, %d\n", visited[pos], matrix[u * dim + v], dist[n*dim+u], dist[pos], v , u);
                if (!visited[pos] && matrix[u * dim + v] && dist[n * dim + u] != dim &&
                    dist[n * dim + u] + matrix[u * dim + v] < dist[pos]) {
                    dist[pos] = dist[n * dim + u] + matrix[u * dim + v];
                }
            }
        }

        // print
        
        printf("Vertiex\t Distance\n");
        for (int i = 0; i < dim; i++)
        {
            printf("Vertiex: %d Distance: %d\n", i, dist[n * dim + i]);
        }
        
       free(visited);
    }

    return dist;
}


// Reads the file and allocates it to memory
int* readFile(FILE *fp, int *dim) {
    int num;
    int nelements = 0;
    int *ptr;
    int *matrix;

    // Reads the dimensions
    fread(&num, sizeof(int), 1, fp);
    *dim = num;

    // initialize the matrix
    matrix = initMatrix(*dim);

    // Reads the elements
    while (fread(&num, sizeof(int), 1, fp) != 0) {
        matrix[nelements++] = num;
    }
    
    return matrix;
}

// Reads the file and allocates it to memory
int readFileDims(FILE *fp) {
    int num;


    // Reads the dimensions
    fread(&num, sizeof(int), 1, fp);
    return num;
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

// Prints the distance matrix
void printDistance(int *root_dist, int dim, int nelements) {
    printf("Distance matrix: \n");
    for (int i = 0; i < nelements; i++)
    {
        printf("%d ", root_dist[i]);
        if (i % dim == dim-1) printf("\n");
    }
}

// Checks that the file can be opened
void fileCheck(FILE *fp) {
    if (fp == NULL) {
        perror("Error: This file cannot be accessed.");
        exit(EXIT_FAILURE);
    }
    return;
}

void memory_check(int *matrix, char *msg) {
    if (matrix == NULL) {
        
        fprintf(stderr, "Memory allocation error\n");
        printf("Error with %s\n", msg);
        exit(EXIT_FAILURE);
    }
}

void mpi_error_check(int mpierror){
    if(mpierror != MPI_SUCCESS){
        printf("Error mpi failure\n");
        fprintf(stderr, "Error with mpi\n");
        exit(EXIT_FAILURE);
    }
}

/** MPI All reduce
 *  Accessing a reduced result across all processes, similar to a mpi_reduce and mpi_broadcast
 *  @param local_min value that each process wants to reduce
 *  @param global_min global value that processes will receive
 *  @param 2 number of elements
 *  @param MPI data type
 *  @param op operation to be performed for reduction
 *  @param communicator
 */
/*
mpierror = MPI_Allreduce(local_min, global_min, 2, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
*/

/** MPI Gather
 *  Collects all the results from the processes into the root process.
 */
//mpierror = MPI_Gather(local_matrix, local_n, MPI_INT, root_matrix, local_n, MPI_INT, root, MPI_COMM_WORLD);