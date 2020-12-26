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
    int world_size = 0;

    if (argc != 2)
    {

        printf("Usage: %s <N>\n", argv[0]);

        return 0;
    }
    srand(time(NULL));
    size_t n = atoi(argv[1]);
    size_t nOverK = n /* / world_size */;

    double *a = allocarray1D(n * n);
    double *b = allocarray1D(n);
    double *x = allocarray1D(n);
    double *x_partial = allocarray1D(nOverK);
    double *xseq = allocarray1D(n);

    double *a_partial = allocarray1D(n * nOverK);

    if (a == NULL || b == NULL || x == NULL || xseq == NULL || x_partial == NULL)
    {
        printf("Allocation failed\n");
        return 0;
    }
    // Process 0 creates A matrix.
    // if (taskid == 0)
    // {
    //     fullfillArrayWithRandomNumbers(a, n * n);
    //     // Process 0 produces the b
    //     fullfillArrayWithRandomNumbers(b, n);
    // }

    // Process 0 sends a_partial to everyone
    // if (!(world_size == 1 && n == 64000))
    // {
    //     MPI_Scatter(a, n * nOverK, MPI_DOUBLE, a_partial, n * nOverK, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    // }

    SequentialMatrixMultiply(n, a, b, xseq);
    // check difference between x and xseq using OpenMP
    //print_1d_arr(exec_times, world_size);
    // print_1d_arr(xseq, n);
    double max_exec, min_exec, avg_exec;
    min_exec = 1000;

    avg_exec = avg_exec / world_size;

    double time_end_openmp, openmp_exec_time, min_exec_time, max_exec_time, avg_exec_time;
    max_exec_time = 0;
    max_exec_time = 1000;
    long double l2_norm = 0;
    size_t numberOfThreads = 0;
    size_t r = 0;
    double *diff_vector = allocarray1D(n);
    size_t nrepeat = 100000;
    double time_start_openmp = omp_get_wtime();
    if (world_size == 1)
    {
        printf("%d times repating\n", nrepeat);
#pragma omp parallel reduction(+ \
                               : l2_norm)
        {
            numberOfThreads = omp_get_num_threads();
            for (int r = 0; r < nrepeat; r++)
            {
                l2_norm = 0;
#pragma omp for
                for (int i = 0; i < n; i++)
                {
                    double local_diff = x[i] - xseq[i];
                    diff_vector[i] = local_diff;
                    l2_norm += (local_diff * local_diff);
                }
            }
        }
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
        printf("OPENMP: %d %ld %f %.12e\n", n, numberOfThreads, openmp_exec_time, openmp_exec_time, l2_norm);
    }
    printf("MIN_AVG_MAX: %d %d %f %f %f\n", n, world_size, min_exec, max_exec, avg_exec);
    printf("MPI: %d %d %f %.12Lf %.12e\n", n, world_size, max_exec, l2_norm, l2_norm);
    totalMemUsage = totalMemUsage / (1024 * 1024 * 1024);
    printf("TOTALMEMUSAGE: %zu\n", totalMemUsage);

    //printf("process: %d %d %d %f %.12e\n", taskid, n, world_size, parallel_exec_time, l2_norm);
    //printf("%d %ld %f %.12e\n", n, numberOfThreads, openmp_exec_time, l2_norm);
    return 0;
}
