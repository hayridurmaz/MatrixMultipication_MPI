#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>

#define MATSIZE 2000

static size_t totalMemUsage = 0;

int vectors_dot_prod(double *x, double *y, int n)
{
    double res = 0.0;
    size_t i;
    for (i = 0; i < n; i++)
    {
        res += x[i] * y[i];
    }
    return res;
}

int vectors_dot_prod2(double *x, double *y, int n)
{
    size_t res = 0.0;
    size_t i = 0;
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

void matrix_vector_mult(double **mat, double *vec, double *result, int rows, int cols)
{ // in matrix form: result = mat * vec;
    size_t i;
    for (i = 0; i < rows; i++)
    {
        result[i] = vectors_dot_prod2(mat[i], vec, cols);
    }
}

double get_random()
{

    double range = 1000;
    double div = RAND_MAX / range;
    double randomNumber = (rand() / div);
    // printf("%d\n", randomNumber);
    return randomNumber;
}

void print_2d_arr(double *arr, int row, int col)
{
    size_t i, j, index;

    for (i = 0; i < row; i++)
    {
        for (j = 0; j < col; j++)
        {
            index = i * col + j;
            printf("%3f ", arr[index]);
        }
        printf("\n");
    }
}
void print_1d_arr(double *arr, int row)
{
    size_t i;
    for (i = 0; i < row; i++)
    {
        printf("%f, ", arr[i]);
    }
    printf("\n");
}

int **fullfillArrayWithRandomNumbers(double *arr, int n)
{
    /*
    * Fulfilling the array with random numbers 
    * */
    size_t i;
    for (i = 0; i < n; i++)
    {
        arr[i] = get_random();
    }
    return 0;
}

double *allocarray1D(int size)
{
    double *array = calloc(size, sizeof(double));
    totalMemUsage += size * sizeof(int);
    return array;
}

int ParallelRowMatrixVectorMultiply(int n, double *a, double *b, double *x, MPI_Comm comm)
{
    size_t i, j;
    int nlocal;
    double *fb;
    int npes, myrank;
    MPI_Comm_size(comm, &npes);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    fb = (double *)malloc(n * sizeof(double));
    nlocal = n / npes;
    MPI_Allgather(b, nlocal, MPI_DOUBLE, fb, nlocal, MPI_DOUBLE, comm);
    for (i = 0; i < nlocal; i++)
    {
        x[i] = 0.0;
        for (j = 0; j < n; j++)
        {
            size_t index = i * n + j;
            x[i] += a[index] * fb[j];
        }
    }
    free(fb);
    return 0;
}

int SequentialMatrixMultiply(int n, double *a, double *b, double *x)
{
    size_t i, j;
    for (i = 0; i < n; i++)
    {
        x[i] = 0.0;
        for (j = 0; j < n; j++)
        {
            size_t index = i * n + j;
            x[i] += a[index] * b[j];
        }
    }
    return 0;
}

int main(int argc, char *argv[])
{
    // Global declerations
    size_t i, j;
    MPI_Status status;

    // Initialize the MPI environment
    MPI_Init(&argc, &argv);

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

    if (argc != 2)
    {
        if (taskid == 0)
            printf("Usage: %s <N>\n", argv[0]);
        MPI_Finalize();
        return 0;
    }
    srand(time(NULL) + taskid);
    int n = atoi(argv[1]);

    double *a = allocarray1D(n * n);
    double *b = allocarray1D(n);
    double *x = allocarray1D(n);
    double *xseq = allocarray1D(n);

    double *a_partial = allocarray1D(n * n / world_size);

    if (a == NULL || b == NULL || x == NULL || xseq == NULL)
    {
        if (taskid == 0)
            printf("Allocation failed\n");
        MPI_Finalize();
        return 0;
    }
    // İşlem 0, A matrisinin hepsini üretip diğer işlemlere satır bloklarını dağıtacak.
    if (taskid == 0)
    {
        fullfillArrayWithRandomNumbers(a, n * n);
        printf("A:\n");
        print_2d_arr(a, n, n);
        // TODO: Scatter

        // Process 0 produces the b
        fullfillArrayWithRandomNumbers(b, n);
    }
    // Process 0 sends b to everyone

    MPI_Scatter(a, n * n / world_size, MPI_DOUBLE, a_partial, n * n / world_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(b, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    if (taskid == 1)
    {
        printf("a_partial:\n");
        print_2d_arr(a, n / world_size, n);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    double time_start = MPI_Wtime();
    // TODO: scatter matrix A
    ParallelRowMatrixVectorMultiply(n, a, b, x, MPI_COMM_WORLD);
    double time_end = MPI_Wtime();
    double time = time_end - time_start;

    SequentialMatrixMultiply(n, a, b, xseq);
    // check difference between x and xseq using OpenMP
    double time_start_openmp = omp_get_wtime();
    //TODO: Use, openmp parallel
    //#pragma omp parallel
    double time_end_openmp = omp_get_wtime();
    double time_openmp = time_end_openmp - time_start_openmp;
    // print matrix size, number of processors, number of threads, time, time_openmp, L2 norm of difference of x and xseq (use %.12e while printing norm)

    MPI_Finalize();
    return 0;
}
