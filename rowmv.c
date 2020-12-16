#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include <math.h>

#define MATSIZE 2000

static size_t totalMemUsage = 0;

size_t vectors_dot_prod(double *x, double *y, size_t n)
{
    double res = 0.0;
    size_t i;
    for (i = 0; i < n; i++)
    {
        res += x[i] * y[i];
    }
    return res;
}

size_t vectors_dot_prod2(double *x, double *y, size_t n)
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

void matrix_vector_mult(double **mat, double *vec, double *result, size_t rows, size_t cols)
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

void print_2d_arr(double *arr, size_t row, size_t col)
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
void print_1d_arr(double *arr, size_t row)
{
    size_t i;
    for (i = 0; i < row; i++)
    {
        printf("%f, ", arr[i]);
    }
    printf("\n");
}

size_t **fullfillArrayWithRandomNumbers(double *arr, size_t n)
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

double *allocarray1D(size_t size)
{
    double *array = calloc(size, sizeof(double));
    totalMemUsage = totalMemUsage + size * sizeof(double);
    return array;
}

size_t ParallelRowMatrixVectorMultiply(size_t n, double *a, double *b, double *x, MPI_Comm comm)
{
    size_t i, j;
    size_t nlocal;
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

size_t ParallelRowMatrixVectorMultiply_WithoutAllgather(size_t n, double *a, double *b, double *x_partial, double *x, MPI_Comm comm)
{

    // Process 0 sends b to everyone
    MPI_Bcast(b, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    size_t i, j;
    size_t nlocal;
    // double *fb;
    int npes, myrank;
    MPI_Comm_size(comm, &npes);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    // fb = (double *)malloc(n * sizeof(double));
    nlocal = n / npes;
    // MPI_Allgather(b, nlocal, MPI_DOUBLE, fb, nlocal, MPI_DOUBLE, comm);
    for (i = 0; i < nlocal; i++)
    {
        x_partial[i] = 0.0;
        for (j = 0; j < n; j++)
        {
            size_t index = i * n + j;
            // printf("%f x %f\n", a[index], b[j]);
            x_partial[i] += a[index] * b[j];
        }
    }
    // free(b);

    // Process 0 gathers x_partials to create x
    MPI_Gather(x_partial, nlocal, MPI_DOUBLE, x, nlocal, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    return 0;
}

size_t SequentialMatrixMultiply(size_t n, double *a, double *b, double *x)
{
    size_t i, j;
    for (i = 0; i < n; i++)
    {
        x[i] = 0.0;
        for (j = 0; j < n; j++)
        {
            size_t index = i * n + j;
            // printf("%f x %f\n", a[index], b[j]);
            x[i] += a[index] * b[j];
        }
    }
    return 0;
}

int main(int argc, char *argv[])
{
    // Global declerations
    size_t i;
    // MPI_Status status;

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
    size_t n = atoi(argv[1]);
    size_t nOverK = n / world_size;

    double *a = allocarray1D(n * n);
    double *b = allocarray1D(n);
    double *x = allocarray1D(n);
    double *x_partial = allocarray1D(nOverK);
    double *xseq = allocarray1D(n);

    double *a_partial = allocarray1D(n * nOverK);

    if (a == NULL || b == NULL || x == NULL || xseq == NULL || x_partial == NULL)
    {
        if (taskid == 0)
            printf("Allocation failed\n");
        MPI_Finalize();
        return 0;
    }
    // Process 0 creates A matrix.
    if (taskid == 0)
    {
        fullfillArrayWithRandomNumbers(a, n * n);
        // Process 0 produces the b
        fullfillArrayWithRandomNumbers(b, n);
    }

    // Process 0 sends a_partial to everyone
    if (!(world_size == 1 && n == 64000))
    {
        MPI_Scatter(a, n * nOverK, MPI_DOUBLE, a_partial, n * nOverK, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    double time_start = MPI_Wtime();
    ParallelRowMatrixVectorMultiply_WithoutAllgather(n, a_partial, b, x_partial, x, MPI_COMM_WORLD);
    double time_end = MPI_Wtime();
    double parallel_exec_time = time_end - time_start;

    double *exec_times = allocarray1D(world_size);
    // Process 0 gathers x_partials to create x
    MPI_Gather(&parallel_exec_time, 1, MPI_DOUBLE, exec_times, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    // print_1d_arr(x, n);

    if (taskid == 0)
    {
        SequentialMatrixMultiply(n, a, b, xseq);
        // check difference between x and xseq using OpenMP
        //print_1d_arr(exec_times, world_size);
        // print_1d_arr(xseq, n);
        double max_exec, min_exec, avg_exec;
        min_exec = 1000;
        for (i = 0; i < world_size; i++)
        {
            if (max_exec < exec_times[i])
            {
                max_exec = exec_times[i];
            }
            if (min_exec > exec_times[i])
            {
                min_exec = exec_times[i];
            }
            avg_exec += exec_times[i];
        }
        avg_exec = avg_exec / world_size;

        double time_start_openmp = omp_get_wtime();
        double time_end_openmp, openmp_exec_time, min_exec_time, max_exec_time, avg_exec_time;
        max_exec_time = 0;
        max_exec_time = 1000;
        long double l2_norm = 0;
        size_t numberOfThreads = 0;
        size_t r = 0;
        double *diff_vector = allocarray1D(n);
        size_t nrepeat = 100000;

        if (world_size == 1)
        {
            printf("%d times repating\n", nrepeat);
            // for (r = 0; r < nrepeat; r++)
            // {
            l2_norm = 0;
#pragma omp parallel for private(i)
            for (i = 0; i < n; i++)
            {
                numberOfThreads = omp_get_num_threads();
                double local_diff = x[i] - xseq[i];
                diff_vector[i] = local_diff;
                l2_norm += (local_diff * local_diff);
            }
            //}
        }
        else
        {
#pragma omp parallel
            {
                numberOfThreads = omp_get_num_threads();
#pragma omp parallel for private(i)
                for (i = 0; i < n; i++)
                {
                    double local_diff = x[i] - xseq[i];
                    diff_vector[i] = local_diff;
                    l2_norm += (local_diff * local_diff);
                }
            }
        }
        l2_norm = sqrt(l2_norm);
        time_end_openmp = omp_get_wtime();
        openmp_exec_time = time_end_openmp - time_start_openmp;
        // print matrix size, number of processors, number of threads, time, time_openmp, L2 norm of difference of x and xseq (use %.12e while printing norm)
        if (world_size == 1)
        {
            printf("OPENMP: %d %ld %f %.12e\n", n, numberOfThreads, openmp_exec_time, l2_norm);
        }
        printf("MIN_AVG_MAX: %d %d %f %f %f\n", n, world_size, min_exec, max_exec, avg_exec);
        printf("MPI: %d %d %f %.12e\n", n, world_size, max_exec, l2_norm);
        totalMemUsage = totalMemUsage / (1024 * 1024 * 1024);
        printf("TOTALMEMUSAGE: %zu\n", totalMemUsage);

        //printf("process: %d %d %d %f %.12e\n", taskid, n, world_size, parallel_exec_time, l2_norm);
        //printf("%d %ld %f %.12e\n", n, numberOfThreads, openmp_exec_time, l2_norm);
    }
    MPI_Finalize();
    return 0;
}
