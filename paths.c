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

    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);


    if (myRank == root) {

        //Conditional checking for invalid command line argument. Could potentially check for greater then 3 or 4 based on the -f flag, however this may be overkill.
        if (argc < 2 )
        {
            printf("An invalid number of arguments has be inputted. \nPlease check your command line arguments.\n");
            exit(0);
        }
        char * file = getFileName( argc, argv );
        printf("The file name is:\t%s\n%d\n", file,argc);

        // Open file
        fp = fopen(file, "rb");
        fileCheck(fp);

        // Reads the file data and sets the matrix array
        root_matrix = readFile(fp, &dim);
        fclose(fp);

        nelements = dim * dim;

        for (int i = 0; i < nelements; i++)
        {
            printf("%d ", root_matrix[i]);
            if (i % 4 == 3) printf("\n");
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
    // broadcasts the number of dimensions
    mpierror = MPI_Bcast(&dim, 1, MPI_INT, root, MPI_COMM_WORLD); mpi_error_check(mpierror);
    printf("DIM: %d %d\n", myRank, dim);

    local_n = (dim * dim) / numProcs; // p / n
    printf("local:%d\n", local_n);

    local_matrix = initMatrix(dim * local_n * sizeof(int));
    local_dist = initMatrix(local_n * sizeof(int));

    //mpierror = MPI_Scatter(root_matrix, local_n, MPI_INT, local_matrix, local_n, MPI_INT, root, MPI_COMM_WORLD);
    //mpi_error_check(mpierror);

    for (int i = 0; i < local_n; i++)
    {
        printf("My process: %d Matrix: %d\n", myRank, local_matrix[i]);
        local_matrix[i]++;
    }

    // Reduce ALL
    //int local_min[2] = { 3, 2 };
    //int global_min[2];


    //mpierror = MPI_Allreduce(local_min, global_min, 2, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    //printf("A Min: %d %d %d %d\n", local_min[0], local_min[1], global_min[0], global_min[1]);
    
    //mpierror = MPI_Gather(local_matrix, local_n, MPI_INT, root_matrix, local_n, MPI_INT, root, MPI_COMM_WORLD);
    //mpi_error_check(mpierror);
    

    if (myRank == root && 0) {
        dist = dijkstra(root_matrix, dim);
        
        for (int i = 0; i < nelements; i++)
        {
            printf("%d ", root_matrix[i]);
            if (i % 4 == 3) printf("\n");
        }
        free(root_matrix);
    }

    //mpierror = MPI_Gather(matrix, partition, MPI_INT, root_matrix, partition, MPI_INT, root, MPI_COMM_WORLD);
    //mpi_error_check(mpierror);

    // Does dijkstra all shortest paths,
    
    //free(root_dist);

    free(local_matrix);
    free(local_dist);
    MPI_Finalize();
    return 0;
}

// Initializes a 2d matrix
int* initMatrix(int dim) {
    int *matrix;

    // Allocate memory depending on dimensions
    matrix = (int *) malloc(sizeof(int) * dim * dim);
    memory_check(matrix);
    return matrix;
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

            for (int v = 0; v < dim; v++)
            {
                pos = n * dim + v;
                //printf("%d, %d, %d, %d, %d, %d\n", visited[pos], matrix[u * dim + v], dist[n*dim+u], dist[pos], v , u);
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
            pos = n * dim + i;
            printf("Vertiex: %d Distance: %d\n", i, dist[pos]);
        }
    }

    // Prints the matrix
    printf("%d\n",dim);
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            pos = i * dim + j;
            printf("%d ",matrix[pos]);
        }
        printf("\n");
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

// Checks that the file can be opened
void fileCheck(FILE *fp) {
    if (fp == NULL) {
        perror("Error: This file cannot be accessed.");
        exit(EXIT_FAILURE);
    }
    return;
}

void memory_check(int *matrix) {
    if (matrix == NULL) {
        printf("Memory allocation error\n");
        fprintf(stderr, "Memory allocation error\n");
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