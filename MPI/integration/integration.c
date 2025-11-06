//
// Created by dominykas on 11/2/25.
//
#include "integration.h"

double func(double x)
{
//    double t = MPI_Wtime(); /* artificially increase the single func evaluation time */
//    Introduce work imballance by sleeping more given larger x
//    while (MPI_Wtime()-t <= IDLETIME*x*x);
    return pow(sin(x),2)/(1+x*x+pow(cos(x),3));
}

double shitty_integrate(double (*f)(double x),
                        double x_start,
                        double x_end,
                        int maxSteps)
{
    int myRank, numProcs;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

    double sum = 0.0;
    double x[2], y[2];

    double stepSize = (x_end - x_start)/(double)maxSteps;
    int step;
    int nextRank = 1;

    double* array = myRank == 0 ? (double*)malloc(sizeof(double) * numProcs) : NULL;

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Status status;

    if (myRank == 0)
    {
        // I am the controller, distribute the work
        for (step = 0; step < maxSteps; step++)
        {
            x[0] = x_start + stepSize*step;
            x[1] = x_start + stepSize*(step+1);
            nextRank = step % (numProcs-1) + 1;
            // Send the work
            MPI_Send(x, 2, MPI_DOUBLE, nextRank, TAG_WORK, MPI_COMM_WORLD);
            // Receive the result
            MPI_Recv(y, 2, MPI_DOUBLE, nextRank, TAG_WORK, MPI_COMM_WORLD,
                     MPI_STATUS_IGNORE);
            sum += stepSize*0.5*(y[0]+y[1]);
        }
        // Signal workers to stop by sending empty messages with tag TAG_END
        for (nextRank = 1; nextRank < numProcs; nextRank++)
            MPI_Send(&nextRank, 0, MPI_INT, nextRank, TAG_END, MPI_COMM_WORLD);
    }
    else
    {
        while (1)
        {
            // I am a worker, wait for work

            // Receive the left and right points of the trapezoid and compute
            // the corresponding function values. If the tag is TAG_END, don't
            // compute but exit.
            MPI_Recv(x, 2, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD,
                     &status);
            if (status.MPI_TAG == TAG_END) break;
            y[0] = f(x[0]);
            y[1] = f(x[1]);
            // Send back the computed result
            MPI_Send(y, 2, MPI_DOUBLE, 0, TAG_WORK, MPI_COMM_WORLD);
        }
    }
    return sum;
}

double integrate(double (*f)(double x),
                 double x_start,
                 double x_end,
                 int maxSteps)
{
    int myRank, numProcs;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

    double sum = 0.0;
    double x[2], y[2];

    double stepSize = (x_end - x_start)/(double)maxSteps;
    int step;
    int nextRank = 1;

    double* array = myRank == 0 ? (double*)malloc(sizeof(double) * numProcs) : NULL;

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Status status;
    double sumA = 0;
    int num_of_steps_local = maxSteps / numProcs;
    int begin = num_of_steps_local * myRank;
    int end = myRank + 1 == numProcs ? maxSteps : num_of_steps_local + begin;
    x_start = begin * stepSize;
    for(int i = begin; i < end; i ++){
       y[0] = f(x_start);
       y[1] = f(x_start + stepSize);
       sumA += stepSize*0.5*(y[0]+y[1]);
       x_start += stepSize;
    }

    //now, broadcast
    MPI_Gather(&sumA, 1, MPI_DOUBLE, array, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if(myRank == 0){
        sumA = 0;
        for(int i = 0; i < numProcs; i ++){
            sumA += array[i];
        }
        free(array);
    }
    return sumA;
}

int integration_dispatch(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    int i = 0;
    int myRank, numProcs;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

    // Integration domain is [0, 1]
    double x0 = 0.0, x1 = 1.0;
    int maxSteps = 100000;

/* interactive reading not possosible when submit to a batch queue with srun */
//    if (myRank == 0)
//    {
//        if (argc > 1)
//        {
//            maxSteps = atoi(argv[1]);
//            if (maxSteps < 1) MPI_Abort(MPI_COMM_WORLD, 1);
//        }
//        printf("Integrating from %lf to %lf in %i steps\n",
//            x0, x1, maxSteps);
//    }

    // Synchronize before making performance measurements
    MPI_Barrier(MPI_COMM_WORLD);

    double startTime = MPI_Wtime();

    double Intgr_val = integrate(func, x0, x1, maxSteps);

    double stopTime = MPI_Wtime();

    if (myRank == 0)
        printf("\nIntegral value = %lf\nComputation took %.7lf seconds\n",
               Intgr_val, stopTime-startTime);

    MPI_Finalize();
    return 0;
}
