//
// Created by dominykas on 11/5/25.
//
/************************************************************
 * Program to solve a finite difference
 * discretization of the screened Poisson equation:
 * (d2/dx2)u + (d2/dy2)u - alpha u = f
 * with zero Dirichlet boundary condition using the iterative
 * Jacobi method with overrelaxation.
 *
 * RHS (source) function
 *   f(x,y) = -alpha*(1-x^2)(1-y^2)-2*[(1-x^2)+(1-y^2)]
 *
 * Analytical solution to the PDE
 *   u(x,y) = (1-x^2)(1-y^2)
 *
 *
 * Input : n     - grid dimension in x direction
 *         m     - grid dimension in y direction
 *         alpha - constant (always greater than 0.0)
 *         tol   - error tolerance for the iterative solver
 *         relax - Successice Overrelaxation parameter
 *         mits  - maximum iterations for the iterative solver
 *
 * On output
 *       : u(n,m)       - Dependent variable (solution)
 *       : f(n,m,alpha) - Right hand side function
 *
 *************************************************************/


#include <f2fs_fs.h>
#include "jacobi.h"
/*************************************************************
 * Performs a halo swap
 *************************************************************/

/**
 * @brief Exchange halo (boundary) rows between neighboring MPI ranks.
 *
 *
 * @param[in]  row_size   Number of elements in each row (width of the domain).
 * @param[in]  local_rows Number of interior rows owned by the current rank
 *                        (excluding halo rows).
 * @param[in,out] src     Local solution array (unique to a process)
 * @param[in]  prevRank   Rank of the previous (upper) neighbor in the MPI
 *                        communicator. Use MPI_PROC_NULL if none.
 * @param[in]  nextRank   Rank of the next (lower) neighbor in the MPI
 *                        communicator. Use MPI_PROC_NULL if none.
 */
void swap_halos(int row_size, int local_rows, double *src, int prevRank, int nextRank) {
    #define ROW(YY) &src[(YY)*row_size]
    // TODO: Send last data row to nextRank
    //       receive first halo data from prevRank
    MPI_Sendrecv(
            ROW(local_rows), row_size, MPI_DOUBLE, nextRank, 0,
            ROW(0),          row_size, MPI_DOUBLE, prevRank, 0,
            MPI_COMM_WORLD, MPI_STATUS_IGNORE
    );

    // Send first interior row to prevRank, receive bottom halo from nextRank
    MPI_Sendrecv(
            ROW(1),              row_size, MPI_DOUBLE, prevRank, 1,
            ROW(local_rows + 1), row_size, MPI_DOUBLE, nextRank, 1,
            MPI_COMM_WORLD, MPI_STATUS_IGNORE
    );
}


/*************************************************************
 * Performs one iteration of the Jacobi method and computes
 * the residual value.
 *
 * NOTE: u(0,*), u(maxXCount-1,*), u(*,0) and u(*,maxYCount-1)
 * are BOUNDARIES and therefore not part of the solution.
 *************************************************************/
double one_jacobi_iteration(double xStart, double yStart, int maxXCount, int maxYCount,double *src, double *dst,
                            double deltaX, double deltaY, double alpha, double omega) {
#define SRC(XX,YY) src[(YY)*maxXCount+(XX)]
#define DST(XX,YY) dst[(YY)*maxXCount+(XX)]
    int x, y;
    double fX, fY;
    double error = 0.0;
    double updateVal;
    double f;
    // Coefficients
    double cx = 1.0/(deltaX*deltaX);
    double cy = 1.0/(deltaY*deltaY);
    double cc = -2.0*cx-2.0*cy-alpha;

    for (y = 1; y < (maxYCount-1); y++)
    {
        fY = yStart + (y-1)*deltaY;
        for (x = 1; x < (maxXCount-1); x++)
        {
            fX = xStart + (x-1)*deltaX;
            f = -alpha*(1.0-fX*fX)*(1.0-fY*fY)-2.0*(1.0-fX*fX)-2.0*(1.0-fY*fY);
            updateVal = ((SRC(x-1,y)+SRC(x+1,y))*cx +
                         (SRC(x,y-1)+SRC(x,y+1))*cy +
                         SRC(x,y)*cc - f)/cc;
            DST(x,y) = SRC(x,y) - omega*updateVal;
            error += updateVal*updateVal;
        }
    }
   return sqrt(error)/((maxXCount-2)*(maxYCount-2));
}

double one_jacobi_iteration_serial(double xStart, double yStart,
                            int maxXCount, int maxYCount,
                            double *src, double *dst,
                            double deltaX, double deltaY,
                            double alpha, double omega)
{
#define SRC(XX,YY) src[(YY)*maxXCount+(XX)]
#define DST(XX,YY) dst[(YY)*maxXCount+(XX)]
    int x, y;
    double fX, fY;
    double error = 0.0;
    double updateVal;
    double f;
    // Coefficients
    double cx = 1.0/(deltaX*deltaX);
    double cy = 1.0/(deltaY*deltaY);
    double cc = -2.0*cx-2.0*cy-alpha;

    for (y = 1; y < (maxYCount-1); y++)
    {
        fY = yStart + (y-1)*deltaY;
        for (x = 1; x < (maxXCount-1); x++)
        {
            fX = xStart + (x-1)*deltaX;
            f = -alpha*(1.0-fX*fX)*(1.0-fY*fY)-2.0*(1.0-fX*fX)-2.0*(1.0-fY*fY);
            updateVal = ((SRC(x-1,y)+SRC(x+1,y))*cx +
                         (SRC(x,y-1)+SRC(x,y+1))*cy +
                         SRC(x,y)*cc - f)/cc;
            DST(x,y) = SRC(x,y) - omega*updateVal;
            error += updateVal*updateVal;
        }
    }
    return sqrt(error)/((maxXCount-2)*(maxYCount-2));
}

/**********************************************************
 * Checks the error between numerical and exact solutions
 **********************************************************/
double checkSolution(double xStart, double yStart,
                     int maxXCount, int maxYCount,
                     double *u,
                     double deltaX, double deltaY,
                     double apha)
{
#define U(XX,YY) u[(YY)*maxXCount+(XX)]
    int x, y;
    double fX, fY;
    double localError, error = 0.0;

    for (y = 1; y < (maxYCount-1); y++)
    {
        fY = yStart + (y-1)*deltaY;
        for (x = 1; x < (maxXCount-1); x++)
        {
            fX = xStart + (x-1)*deltaX;
            localError = U(x,y) - (1.0-fX*fX)*(1.0-fY*fY);
            error += localError*localError;
        }
    }
    return sqrt(error)/((maxXCount-2)*(maxYCount-2));
}

double checkSolutionSerial(double xStart, double yStart,
                     int maxXCount, int maxYCount,
                     double *u,
                     double deltaX, double deltaY,
                     double alpha)
{
#define U(XX,YY) u[(YY)*maxXCount+(XX)]
    int x, y;
    double fX, fY;
    double localError, error = 0.0;

    for (y = 1; y < (maxYCount-1); y++)
    {
        fY = yStart + (y-1)*deltaY;
        for (x = 1; x < (maxXCount-1); x++)
        {
            fX = xStart + (x-1)*deltaX;
            localError = U(x,y) - (1.0-fX*fX)*(1.0-fY*fY);
            error += localError*localError;
        }
    }
    return sqrt(error)/((maxXCount-2)*(maxYCount-2));
}

int jacobi(int argc, char **argv)
{
    int n, m, mits;
    double alpha, tol, relax;
    double maxAcceptableError;
    double localError, error;
    double *u, *u_old, *tmp;
    double *u_big;
    int allocCount;
    int iterationCount, maxIterationCount;

    int myRank, numProcs;
    int prevRank, nextRank;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

    prevRank = (myRank > 0) ? (myRank-1) : MPI_PROC_NULL;
    nextRank = (myRank < numProcs-1) ? (myRank+1) : MPI_PROC_NULL;

    if (myRank == 0)
    {
//        printf("Input n,m - grid dimension in x,y direction:\n");
//        scanf("%d,%d", &n, &m);
//        printf("Input alpha - Helmholtz constant:\n");
//        scanf("%lf", &alpha);
//        printf("Input relax - successive over-relaxation parameter:\n");
//        scanf("%lf", &relax);
//        printf("Input tol - error tolerance for the iterrative solver:\n");
//        scanf("%lf", &tol);
//        printf("Input mits - maximum solver iterations:\n");
//        scanf("%d", &mits);
//        printf("-> %d, %d, %g, %g, %g, %d\n", n, m, alpha, relax, tol, mits);
//        printf("initialized with: %d %d %lf", n, m, tol);
        m = 16000;
        n = 2048;
        alpha = 0.8;
        relax = 1.0f;
        tol = 1e-13;
        mits = 20;
    }



    if (myRank == 0 && (m % numProcs) != 0)
    {
        fprintf(stderr,
                "ERROR: Processor count does not divide problem size\n"
                "m = %d numProcs = %d\n", m, numProcs);
        MPI_Abort(MPI_COMM_WORLD, 1);
        exit(1);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    double startTime;
    if(myRank == 0) startTime = MPI_Wtime();
    // Broadcast configuration data to all processes
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&m, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&alpha, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&relax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&tol, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&mits, 1, MPI_INT, 0, MPI_COMM_WORLD);

    int myStartRow = myRank * (m / numProcs);
    int myRows = m / numProcs;

    allocCount = (n+2)*(myRows+2);
    // Those two calls also zero the boundary elements
    u = (double*)calloc(sizeof(double), allocCount);
    u_old = (double*)calloc(sizeof(double), allocCount);
    if (u == NULL || u_old == NULL)
    {
        fprintf(stderr, "[%d] Not enough memory for two %ix%i matrices\n",
                myRank, n+2, myRows+2);
        MPI_Abort(MPI_COMM_WORLD, 1);
        exit(1);
    }
    // Master process also allocates an array that holds the entire solution
    if (myRank == 0)
    {
        u_big = (double *)calloc(sizeof(double), (n+2)*(m+2));
        if (u_big == NULL)
        {
            fprintf(stderr, "Not enough memory for the whole solution\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
            exit(1);
        }
    }

    maxIterationCount = mits;
    maxAcceptableError = tol;

    // Solve in [-1, 1] x [-1, 1]
    double xLeft = -1.0, xRight = 1.0;
    double yBottom = -1.0, yUp = 1.0;

    double deltaX = (xRight-xLeft)/(n-1);
    double deltaY = (yUp-yBottom)/(m-1);

    double myYBottom = yBottom + myStartRow*deltaY;

    // Distribute the initial solution approximation to all processes
    if (myRank == 0){
        MPI_Scatter(u_big + (n+2), myRows*(n+2), MPI_DOUBLE, u_old + (n+2), myRows*(n+2), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
    else{
        MPI_Scatter(NULL, myRows*(n+2), MPI_DOUBLE, u_old + (n+2), myRows*(n+2), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }


    iterationCount = 0;
    error = HUGE_VAL;

    /* Iterate as long as it takes to meet the convergence criterion */
    while (iterationCount < maxIterationCount && error > maxAcceptableError)
    {
        //swap hallos
        swap_halos(n+2, myRows, u_old, prevRank, nextRank);

        //call jacobi
        localError = one_jacobi_iteration(xLeft, myYBottom, n+2, myRows+2, u_old, u, deltaX, deltaY, alpha, relax);

        // TODO: Perform global reduction and compute global error value
        error = 0;
        MPI_Reduce(&localError, &error, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Bcast(&error, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        if (myRank == 0)printf("\tError %g\n", error);
        iterationCount++;
        // Swap the buffers
        tmp = u_old;
        u_old = u;
        u = tmp;
    }

    if (myRank == 0)printf("Residual %g\n", error);
    // u_old holds the solution after the most recent buffers swap
    double absoluteError = 0.0;
    localError = checkSolution(xLeft, myYBottom, n+2, myRows+2, u_old, deltaX, deltaY, alpha);
    printf("Local err is %f \n", localError);
    // TODO: Perform global reduction and compute absoluteError value
    MPI_Reduce(&localError, &absoluteError, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (myRank == 0)printf("The error of the iterative solution is %g\n", absoluteError);

    // Gather back the whole solution
    if (myRank == 0){
        MPI_Gather(u_old + (n+2), myRows*(n+2), MPI_DOUBLE, u_big + (n+2), myRows*(n+2), MPI_DOUBLE, 0, MPI_COMM_WORLD);

    }
    else{
        MPI_Gather(u_old + (n+2), myRows*(n+2), MPI_DOUBLE,NULL, myRows*(n+2), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    // Check the solution again to compare to the reduction result
    if (myRank == 0)
    {
        absoluteError = checkSolution(xLeft, yBottom, n+2, m+2, u_big, deltaX, deltaY, alpha);
        absoluteError = sqrt(absoluteError)/(n*m);
        printf("The error of the gathered solution is %g\n", absoluteError);
    }
    double endTime;
    if(myRank == 0) endTime = MPI_Wtime();
    if(myRank == 0) printf("%f", endTime - startTime);
    MPI_Finalize();
    return 0;
}
int jacobi_serial(int argc, char **argv)
{
    int n, m, mits;
    double alpha, tol, relax;
    double maxAcceptableError;
    double error;
    double *u, *u_old, *tmp;
    int allocCount;
    int iterationCount, maxIterationCount;

    printf("Input n,m - grid dimension in x,y direction:\n");
    scanf("%d,%d", &n, &m);
    printf("Input alpha - Helmholtz constant:\n");
    scanf("%lf", &alpha);
    printf("Input relax - successive over-relaxation parameter:\n");
    scanf("%lf", &relax);
    printf("Input tol - error tolerance for the iterrative solver:\n");
    scanf("%lf", &tol);
    printf("Input mits - maximum solver iterations:\n");
    scanf("%d", &mits);
    printf("-> %d, %d, %g, %g, %g, %d\n", n, m, alpha, relax, tol, mits);

    allocCount = (n+2)*(m+2);
    // Those two calls also zero the boundary elements
    u = (double*)calloc(sizeof(double), allocCount);
    u_old = (double*)calloc(sizeof(double), allocCount);
    if (u == NULL || u_old == NULL)
    {
        fprintf(stderr, "Not enough memory for two %ix%i matrices\n", n+2, m+2);
        exit(1);
    }
    maxIterationCount = mits;
    maxAcceptableError = tol;

    // Solve in [-1, 1] x [-1, 1]
    double xLeft = -1.0, xRight = 1.0;
    double yBottom = -1.0, yUp = 1.0;

    double deltaX = (xRight-xLeft)/(n-1);
    double deltaY = (yUp-yBottom)/(m-1);

    iterationCount = 0;
    error = HUGE_VAL;
    double startTime = MPI_Wtime();
    /* Iterate as long as it takes to meet the convergence criterion */
    while (iterationCount < maxIterationCount && error > maxAcceptableError)
    {
        printf("Iteration %i\n", iterationCount);
        error = one_jacobi_iteration_serial(xLeft, yBottom,
                                     n+2, m+2,
                                     u_old, u,
                                     deltaX, deltaY,
                                     alpha, relax);
        printf("\tError %g\n", error);
        iterationCount++;
        // Swap the buffers
        tmp = u_old;
        u_old = u;
        u = tmp;
    }
    printf("Residual %g\n",error);

    // u_old holds the solution after the most recent buffers swap
    double absoluteError = checkSolutionSerial(xLeft, yBottom,
                                         n+2, m+2,
                                         u_old,
                                         deltaX, deltaY,
                                         alpha);
    printf("The error of the iterative solution is %g\n", absoluteError);
    double end = MPI_Wtime();
    printf("%f", end - startTime);
    return 0;
}