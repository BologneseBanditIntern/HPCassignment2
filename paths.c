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
    clock_t timer_start, timer_end;

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
    if (myRank == root) printf("The file name is:\t%s\n%d\n", file,argc);

    // Open file
    fp = fopen(file, "rb");
    fileCheck(fp);

    // Reads the file data and sets the matrix array
    root_matrix = readFile(fp, &dim);
    fclose(fp);

    nelements = dim * dim;

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

    local_matrix = initMatrixP(local_n);

    // Start time of operations
    timer_start = clock();

    /** MPI Scatter
     *  One process divides an array into pieces which are distributed among the processors
     *  The root_matrix is split 
     *  @param local_n number of elements being sent/received per process
    */ 
    mpierror = MPI_Scatter(root_matrix, local_n, MPI_INT, local_matrix, local_n, MPI_INT, root, MPI_COMM_WORLD);
    mpi_error_check(mpierror);

    //printf("Scatter: Rank: %d, Array: %d %d %d %d\n", myRank, local_matrix[0], local_matrix[1], local_matrix[2], local_matrix[3]);

    // Reduce ALL
    //int local_min[3] = { 3, 2 };
    //int global_min[3];

    //mpierror = MPI_Allreduce(local_min, global_min, 2, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    //printf("A Min: %d %d %d %d\n", local_min[0], local_min[1], global_min[0], global_min[1]);
    
    //mpierror = MPI_Gather(local_matrix, local_n, MPI_INT, root_matrix, local_n, MPI_INT, root, MPI_COMM_WORLD);
    //mpi_error_check(mpierror);

    if (myRank == root) {
        //root_dist = dijkstra(root_matrix, dim);
    }

    local_dist = dijkstraP(local_matrix, dim, local_n, myRank, root_matrix);

    /** MPI Gather
     *  Collects all the results from the processes into the root process.
     *  Processors send their elements to the root array to be collected.
     *  The elements are ordered by the rank of the process
     *  local_n is the number of elements received per process
     *  local_dist is the elements stored on the process which collected to the root_dist
     */ 
    mpierror = MPI_Gather(local_dist, local_n, MPI_INT, root_dist, local_n, MPI_INT, root, MPI_COMM_WORLD);
    mpi_error_check(mpierror);

    
    if (myRank == root) {
        // Prints the distance matrix
        printDistance(root_dist, dim, nelements);
    }
    

    // Prints the time of the execution
    timer_end = clock();
    double timespent = (double) (timer_end - timer_start) / CLOCKS_PER_SEC;
    printf("%lf\n", timespent);
    
    free(root_dist);
    free(root_matrix);
    free(dist);
    free(local_matrix);
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

int* dijkstraP(int *matrix, int dim, int local_n, int myRank, int *root_matrix) {

    // output array, holds the shortest distances
    int *dist = initMatrixP(local_n);

    // 1 if shortest distance has been found
    int *visited = initMatrixP(local_n);

    int pos;

    // Initialize all distance values as max (dim) and visited
    for (int i = 0; i < local_n; i++)
    {
        dist[i] = dim;
        visited[i] = 0;
    }
    //printf("My Rank: %d:  %d, %d, %d, %d\n", myRank, matrix[0], matrix[1], matrix[2] ,matrix[3]);
    // iterates through all vertices
    for (int n = 0; n < local_n / dim; n++) {

        // distance from self is 0
        dist[(n * dim + n) + myRank * local_n / dim] = 0;

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
            

            //printf("U: %d\n", u);

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
            for (int j = 0; j < dim; j++)
            {
                printf("%d ", dist[i * dim + j]);
            }
            printf("\n");
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