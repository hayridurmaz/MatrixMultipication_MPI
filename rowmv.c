#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define MATSIZE 8

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

int get_random()
{

    int randomNumber = rand() % 100;
    // printf("%d\n", randomNumber);
    return randomNumber;
}

void print_arr(int **arr, int row, int col)
{
    int i, j;
    for (i = 0; i < row; i++)
    {
        for (j = 0; j < col; j++)
        {
            printf("%3d ", arr[i][j]);
        }
        printf("\n");
    }
}

void print_arr_1d(int *arr, int row)
{
    int i, j;
    for (i = 0; i < row; i++)
    {
        printf("%d, ", arr[i]);
    }
    printf("\n");
}

int **allocarray(int row, int col)
{
    int **array = malloc(row * sizeof(int *));
    for (int i = 0; i < row; i++)
    {
        array[i] = malloc(col * sizeof(int));
        // printf("col=%d - Size of array1[%d]=%d (%d)\n", col, i, sizeof(array1[i]), ((int *)(col * sizeof(int))));
    }
    return array;
}

int **fullfillArrayWithRandomNumbers(int **arr, int row, int col)
{
    /*
    * Fulfilling the array with random numbers 
    * */
    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < col; j++)
        {
            arr[i][j] = get_random();
        }
    }
    return arr;
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

    //General declerations
    int i, j;
    int N = MATSIZE;
    int K, nOverK;
    MPI_Status status;

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

    // Proccess depended inits
    K = world_size;
    nOverK = N / K;
    srand(time(NULL) + taskid);

    // Every proccess creates its part of matrix.
    int **partial_matrix = allocarray(nOverK, N);

    partial_matrix = fullfillArrayWithRandomNumbers(partial_matrix, nOverK, N);
    print_arr(partial_matrix, nOverK, N);

    int **totalArr = NULL;
    if (taskid == 0)
    {
        totalArr = allocarray(N, N);
    }
    // MPI_Gather((arr), N * nOverK, MPI_INT, (totalArr), N * N, MPI_INT, 0,
    //            MPI_COMM_WORLD);
    // free(arr[0]);
    // free(arr);
    if (taskid == 0)
    {
        print_arr(totalArr, N, N);
        // MPI_Send(&(arr[0][0]), N * nOverK, MPI_INT, 1, 111, MPI_COMM_WORLD);
    }
    else if (taskid == 1)
    {
        // int **A = allocarray(nOverK, N);
        // MPI_Recv(&(A[0][0]), nOverK * N, MPI_INT, 0, 111, MPI_COMM_WORLD, &status);
        // if (status.MPI_ERROR == MPI_SUCCESS)
        // {
        //     print_arr(A, nOverK, N);
        // }
        // else
        // {
        //     printf("Rank %d encountered problem at line %d of file %s. \
        //                 Source:%d Tag:%d Error code:%d\n",
        //            taskid, __LINE__, __FILE__,
        //            status.MPI_SOURCE,
        //            status.MPI_TAG,
        //            status.MPI_ERROR);
        // }
        // free(A);
    }
    else
    {
        /* code */
    }

    // Finalize the MPI environment.
    MPI_Finalize();
    return 0;
}