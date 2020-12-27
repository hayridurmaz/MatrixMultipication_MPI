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
    double *b = allocarray1D(n * n);
    double *x = allocarray1D(n * n);

    if (a == NULL || b == NULL || x == NULL)
    {
        printf("Allocation failed\n");
        return 0;
    }

    fullfillArrayWithRandomNumbers(a, n * n);
    fullfillArrayWithRandomNumbers(b, n * n);

    double time_end_openmp, openmp_exec_time;
    size_t numberOfThreads = 0;
    size_t r = 0;
    double *diff_vector = allocarray1D(n);
    double time_start_openmp = omp_get_wtime();
    if (world_size == 1)
    {
#pragma omp parallel
        {
            numberOfThreads = omp_get_num_threads();
#pragma omp for
            for (int i = 0; i < n; i++)
            {
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
            }
        }
    }
    time_end_openmp = omp_get_wtime();
    openmp_exec_time = time_end_openmp - time_start_openmp;
    // print matrix size, number of processors, number of threads, time, time_openmp, L2 norm of difference of x and xseq (use %.12e while printing norm)
    // printf("OPENMP: %d %ld %f %.12e\n", n, numberOfThreads, openmp_exec_time, openmp_exec_time, l2_norm);
    // printf("MIN_AVG_MAX: %d %d %f %f %f\n", n, world_size, min_exec, max_exec, avg_exec);
    // printf("MPI: %d %d %f %.12Lf %.12e\n", n, world_size, max_exec, l2_norm, l2_norm);
    totalMemUsage = totalMemUsage / (1024 * 1024 * 1024);
    printf("TOTALMEMUSAGE: %zu\n", totalMemUsage);

    //printf("process: %d %d %d %f %.12e\n", taskid, n, world_size, parallel_exec_time, l2_norm);
    //printf("%d %ld %f %.12e\n", n, numberOfThreads, openmp_exec_time, l2_norm);
    return 0;
}
