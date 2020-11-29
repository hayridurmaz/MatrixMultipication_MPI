#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>

#define MATSIZE 5

double vectors_dot_prod(const double *x, const double *y, int n)
{
    double res = 0.0;
    int i;
    for (i = 0; i < n; i++)
    {
        res += x[i] * y[i];
    }
    return res;
}

double vectors_dot_prod2(const double *x, const double *y, int n)
{
    double res = 0.0;
    int i = 0;
    for (; i <= n - 4; i += 4)
    {
        res += (x[i] * y[i] +
                x[i + 1] * y[i + 1] +
                x[i + 2] * y[i + 2] +
                x[i + 3] * y[i + 3]);
    }
    for (; i < n; i++)
    {
        res += x[i] * y[i];
    }
    return res;
}

void matrix_vector_mult(const double **mat, const double *vec, double *result, int rows, int cols)
{ // in matrix form: result = mat * vec;
    int i;
    for (i = 0; i < rows; i++)
    {
        result[i] = vectors_dot_prod2(mat[i], vec, cols);
    }
}

double get_random()
{
    double randomNumber = rand() % 100;
    // printf("%f\n", randomNumber);
    return randomNumber;
}

void print_arr(double **arr, int row, int col)
{
    int i, j;
    for (i = 0; i < row; i++)
    {
        for (j = 0; j < col; j++)
        {
            printf("%f ", arr[i][j]);
        }
        printf("\n");
    }
}

void print_arr_1d(int *arr, int row)
{
    int i, j;
    for (i = 0; i < row; i++)
    {
        printf("%f, ", arr[i]);
    }
}

int main(int argc, char *argv[])
{
    // int numtasks,              /* number of tasks in partition */
    //     taskid,                /* a task identifier */
    //     numworkers,            /* number of worker tasks */
    //     source,                /* task id of message source */
    //     dest,                  /* task id of message destination */
    //     mtype,                 /* message type */
    //     rows,                  /* rows of matrix A sent to each worker */
    //     averow, extra, offset, /* used to determine rows sent to each worker */
    //     i, j, k, rc;           /* misc */
    // int tag = 4;
    // MPI_Status status;
    // int N = MATSIZE;
    // double *A,
    //     b[N][1],
    //     x[N][1];

    // printf("Initializing arrays...\n");
    // for (i = 0; i < N; i++)
    //     for (j = 0; j < N; j++)
    //         A[i][j] = get_random();
    // print_arr(A, N, N);

    // MPI_Init(&argc, &argv);
    // MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
    // MPI_Comm_size(MPI_COMM_WORLD, &numtasks);

    //General declerations
    int i;
    int N = MATSIZE;

    // Initialize the MPI environment
    MPI_Init(NULL, NULL);

    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int taskid;
    MPI_Comm_rank(MPI_COMM_WORLD, &taskid);

    // Get the name of the processor
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    // Print off a hello world message
    // printf("Hello world from processor %s, rank %d out of %d processors\n",
    //        processor_name, taskid, world_size);

    if (taskid == 0)
    {
        int *A = calloc(N, sizeof(int));

        for (i = 0; i < N; i++)
            A[i] = get_random();
        MPI_Send(A, N, MPI_INT, 1, 100, MPI_COMM_WORLD);
        free(A);
    }
    else
    {
        int *receivebuffer = calloc(N, sizeof(int));
        MPI_Status status;
        MPI_Recv(receivebuffer, N, MPI_INT, 0, 100, MPI_COMM_WORLD, &status);
        if (status.MPI_ERROR == MPI_SUCCESS)
        {
            printf("receivebuffer:%d\n", &receivebuffer);
            print_arr_1d(receivebuffer, N);
        }
        else
        {
            printf("Rank %d encountered problem at line %d of file %s. \
                    Source:%d Tag:%d Error code:%d\n",
                   taskid, __LINE__, __FILE__,
                   status.MPI_SOURCE,
                   status.MPI_TAG,
                   status.MPI_ERROR);
        }
        free(receivebuffer);
    }

    // Finalize the MPI environment.
    MPI_Finalize();

    // if (taskid == 0)
    // {

    //     // for (int i = 0; i < N; i++)
    //     // {
    //     //     A[i] = (double *)malloc(N * sizeof(double));
    //     // }
    //     A = (double *)malloc(N * sizeof(double));

    //     printf("Initializing arrays...\n");
    //     for (i = 0; i < N; i++)
    //         // for (j = 0; j < N; j++)
    //         A[i] = get_random();
    //     MPI_Send(A, N, MPI_DOUBLE, 1, tag, MPI_COMM_WORLD);
    // }
    // else if (taskid == 1)
    // {
    //     MPI_Recv(A, N, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD,
    //              MPI_STATUS_IGNORE);
    //     printf("I am %d\n", taskid);
    //     print_arr_1d(A, N);
    // }
    // else if (taskid == 2)
    // {
    // }
    // else
    // {
    // }

    // int K = numtasks;

    // if (taskid == 0)
    // {
    // }

    // int *datain, *dataout;

    // datain = (int *)malloc(N * numtasks * sizeof(int));
    // dataout = (int *)malloc(N * numtasks * sizeof(int));

    // MPI_Bcast(A, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // MPI_Allgather(&sa, 2, MPI_INT, &sa, 3, MPI_INT, MPI_COMM_WORLD);

    // for (i = 0; i < N; i++)
    //     b[taskid][i] = get_random();

    // for (i = 0; i < N; i++)
    //     A[taskid][i] = get_random();

    /**************************** master task ************************************/
    // if (taskid == MASTER)
    // {
    //     printf("mpi_mm has started with %d tasks.\n", numtasks);
    //     // printf("Initializing arrays...\n");
    //     for (i = 0; i < NRA; i++)
    //         for (j = 0; j < NCA; j++)
    //             a[i][j] = i + j;
    //     for (i = 0; i < NCA; i++)
    //         for (j = 0; j < NCB; j++)
    //             b[i][j] = i * j;

    //     /* Measure start time */
    //     double start = MPI_Wtime();

    //     /* Send matrix data to the worker tasks */
    //     averow = NRA / numworkers;
    //     extra = NRA % numworkers;
    //     offset = 0;
    //     mtype = FROM_MASTER;
    //     for (dest = 1; dest <= numworkers; dest++)
    //     {
    //         rows = (dest <= extra) ? averow + 1 : averow;
    //         // printf("Sending %d rows to task %d offset=%d\n",rows,dest,offset);
    //         MPI_Send(&offset, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
    //         MPI_Send(&rows, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
    //         MPI_Send(&a[offset][0], rows * NCA, MPI_DOUBLE, dest, mtype,
    //                  MPI_COMM_WORLD);
    //         MPI_Send(&b, NCA * NCB, MPI_DOUBLE, dest, mtype, MPI_COMM_WORLD);
    //         offset = offset + rows;
    //     }

    //     /* Receive results from worker tasks */
    //     mtype = FROM_WORKER;
    //     for (i = 1; i <= numworkers; i++)
    //     {
    //         source = i;
    //         MPI_Recv(&offset, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
    //         MPI_Recv(&rows, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
    //         MPI_Recv(&c[offset][0], rows * NCB, MPI_DOUBLE, source, mtype,
    //                  MPI_COMM_WORLD, &status);
    //         // printf("Received results from task %d\n",source);
    //     }

    //     /* Print results */

    //     printf("******************************************************\n");
    //     printf("Result Matrix:\n");
    //     for (i = 0; i < NRA; i++)
    //     {
    //         printf("\n");
    //         for (j = 0; j < NCB; j++)
    //             printf("%6.2f   ", c[i][j]);
    //     }
    //     printf("\n******************************************************\n");

    //     /* Measure finish time */
    //     double finish = MPI_Wtime();
    //     printf("Done in %f seconds.\n", finish - start);
    // }

    /**************************** worker task ************************************/
    /*     if (taskid > MASTER)
    {
        mtype = FROM_MASTER;
        MPI_Recv(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
        MPI_Recv(&rows, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
        MPI_Recv(&a, rows * NCA, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD, &status);
        MPI_Recv(&b, NCA * NCB, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD, &status);

        for (k = 0; k < NCB; k++)
            for (i = 0; i < rows; i++)
            {
                c[i][k] = 0.0;
                for (j = 0; j < NCA; j++)
                    c[i][k] = c[i][k] + a[i][j] * b[j][k];
            }
        mtype = FROM_WORKER;
        MPI_Send(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
        MPI_Send(&rows, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
        MPI_Send(&c, rows * NCB, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD);
    } */

    // MPI_Finalize();
}