
#include "pingPong-a.h"
// Maximum array size 2^10 = 1024 elements


void ping_pong_b(int argc, char **argv)
{
    // Variables for the process rank and number of processes
    int myRank, numProcs;

    // Initialize MPI
    MPI_Init(&argc, &argv);
    // Find out MPI communicator size and process rank
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Status status;
    // PART B
    srandom(MPI_Wtime()*100000 + myRank*137);
    int numberOfElementsToSend = random() % 100;
    // Allocate an array big enough to hold event the largest message
    int *myArray = (int *)malloc(sizeof(int)*MAX_ARRAY_SIZE);
    int array[256];
    if (myArray == NULL)
    {
        printf("Not enough memory\n");
        exit(1);
    }
    int numberOfElementsReceived;

    // Have only the first process execute the following code
    if (myRank == 0)
    {
        printf("Sending %i elements\n", numberOfElementsToSend);
        // TODO: Send "numberOfElementsToSend" elements
        MPI_Send(array, numberOfElementsToSend, MPI_INT, 1, 0, MPI_COMM_WORLD);
        // TODO: Receive elements
        MPI_Recv(myArray, MAX_ARRAY_SIZE, MPI_INT, 1, 0, MPI_COMM_WORLD, &status);
        // TODO: Store number of elements received in numberOfElementsReceived
        MPI_Get_count(&status, MPI_INT, &numberOfElementsReceived);
    }
    else // myRank == 1
    {
        // TODO: Receive elements
        MPI_Recv(myArray, MAX_ARRAY_SIZE, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
        // TODO: Store number of elements received in numberOfElementsReceived
        MPI_Get_count(&status, MPI_INT, &numberOfElementsReceived);
        printf("Received %i elements\n", numberOfElementsReceived);

        printf("Sending back %i elements\n", numberOfElementsToSend);
        // TODO: Send "numberOfElementsToSend" elements
        MPI_Send(myArray, numberOfElementsToSend, MPI_INT, 0, 0, MPI_COMM_WORLD);

    }
    free(myArray);
    // Finalize MPI
    MPI_Finalize();
}
