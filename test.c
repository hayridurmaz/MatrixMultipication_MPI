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
    // MPI_Status status;

    // Initialize the MPI environment

    if (argc != 2)
    {
        printf("Usage: %s <N>\n", argv[0]);
    }
    srand(time(NULL) + 1);
    size_t n = atoi(argv[1]);

    double *a = allocarray1D(n * n);
    double *b = allocarray1D(n);
    double *x = allocarray1D(n);
    double *xseq = allocarray1D(n);

    fullfillArrayWithRandomNumbers(a, n * n);

    SequentialMatrixMultiply(n, a, b, xseq);

    double time_end_openmp, openmp_exec_time, min_exec_time, max_exec_time, avg_exec_time;
    max_exec_time = 0;
    max_exec_time = 1000;
    long double l2_norm = 0;
    size_t numberOfThreads = 0;
    size_t r = 0;
    double *diff_vector = allocarray1D(n);
    size_t nrepeat = 100000;

    printf("%d times repating\n", nrepeat);
    l2_norm = 0;
    double time_start_openmp = omp_get_wtime();

    printf("%d times repating\n", nrepeat);
#pragma omp parallel
    {
        numberOfThreads = omp_get_num_threads();
        for (int r = 0; r < nrepeat; r++)
        {
            l2_norm = 0;
#pragma omp for reduction(+ \
                          : l2_norm)
            for (int i = 0; i < n; i++)
            {
                double local_diff = x[i] - xseq[i];
                diff_vector[i] = local_diff;
                l2_norm += (local_diff * local_diff);
            }
        }
    }

    time_end_openmp = omp_get_wtime();

    l2_norm = sqrt(l2_norm);

    openmp_exec_time = time_end_openmp - time_start_openmp;
    printf("OPENMP: %d %ld %f %.12e\n", n, numberOfThreads, openmp_exec_time, l2_norm);
    totalMemUsage = totalMemUsage / (1024 * 1024 * 1024);
    printf("TOTALMEMUSAGE: %zu\n", totalMemUsage);
}
