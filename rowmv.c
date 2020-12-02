#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define MATSIZE 2000

static size_t totalMemUsage = 0;

int vectors_dot_prod(int *x, int *y, int n)
{
    int res = 0.0;
    int i;
    for (i = 0; i < n; i++)
    {
        res += x[i] * y[i];
    }
    return res;
}

int vectors_dot_prod2(int *x, int *y, int n)
{
    int res = 0.0;
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

void matrix_vector_mult(int **mat, int *vec, int *result, int rows, int cols)
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

int **allocarray(int row, int col)
{
    int **array = malloc(row * sizeof(int *));
    for (int i = 0; i < row; i++)
    {
        array[i] = malloc(col * sizeof(int));
        totalMemUsage += col * sizeof(int);
        // printf("col=%d - Size of array1[%d]=%d (%d)\n", col, i, sizeof(array1[i]), ((int *)(col * sizeof(int))));
    }
    totalMemUsage += row * sizeof(int *);

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

void print_arr_1d(int *arr, int row)
{
    int i, j;
    for (i = 0; i < row; i++)
    {
        printf("%d, ", arr[i]);
    }
    printf("\n");
}

int *allocarray1D(int row)
{
    int *array = malloc(row * sizeof(int *));
    totalMemUsage += row * sizeof(int *);
    return array;
}

int *fullfillArrayWithRandomNumbers1D(int *arr, int row)
{
    /*
    * Fulfilling the array with random numbers 
    * */
    for (int i = 0; i < row; i++)
    {
        arr[i] = get_random();
    }
    return arr;
}

int main(int argc, char *argv[])
{
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

    // Every proccess creates its part of matrix. (Part 1)
    int **A_partial = allocarray(nOverK, N);

    A_partial = fullfillArrayWithRandomNumbers(A_partial, nOverK, N);
    // print_arr(A_partial, nOverK, N);

    // Every proccess creates its part of B and X vectors. (Part 2)
    int *B_partial = allocarray1D(nOverK);
    B_partial = fullfillArrayWithRandomNumbers1D(B_partial, nOverK);

    int *X_partial = allocarray1D(nOverK);
    X_partial = fullfillArrayWithRandomNumbers1D(X_partial, nOverK);

    //ALLGATHER FOR VECTOR
    int *gatherInput, *B_complete;
    gatherInput = (int *)malloc(N * sizeof(int));
    B_complete = (int *)malloc(N * sizeof(int));
    for (int i = 0; i < N; i++)
        gatherInput[i] = -1;
    for (int i = 0; i < nOverK; i++)
        gatherInput[taskid * nOverK + i] = B_partial[i];

    MPI_Allgather(&(gatherInput[taskid * nOverK]), nOverK, MPI_INT, B_complete, nOverK, MPI_INT, MPI_COMM_WORLD);
    // print_arr(A_partial, nOverK, N);
    // print_arr_1d(B_complete, N);

    MPI_Barrier(MPI_COMM_WORLD);
    double time_start = MPI_Wtime();

    matrix_vector_mult(A_partial, B_complete, X_partial, nOverK, N);
    MPI_Barrier(MPI_COMM_WORLD);
    double time_end = MPI_Wtime();
    double time = time_end - time_start;

    if (taskid == 0)
    {
        totalMemUsage = totalMemUsage / (1024 * 1024);
        printf("%6d %6d %6f %zu\n", N, world_size, time, totalMemUsage);
    }
    // print_arr_1d(X_partial, nOverK);
    // Finalize the MPI environment.
    MPI_Finalize();
    return 0;
}