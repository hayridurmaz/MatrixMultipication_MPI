#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include <math.h>
#include <string.h>

#define min(x, y) (((x) < (y)) ? (x) : (y))

static size_t totalMemUsage = 0;
static size_t numberOfThreads = 0;

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

    double range = 100;
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

int fullfillArrayWithRandomNumbers(double *arr, size_t n)
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

double **allocarray_2D_ArrayOfArrays(int row, int col)
{
    double **array = malloc(row * sizeof(double *));
    size_t i, j;
    for (i = 0; i < row; i++)
    {
        array[i] = malloc(col * sizeof(double));
        totalMemUsage += col * sizeof(double);
        // printf("col=%d - Size of array1[%d]=%d (%d)\n", col, i, sizeof(array1[i]), ((int *)(col * sizeof(int))));
    }
    totalMemUsage += row * sizeof(double *);

    return array;
}

double **fullfillArrayWithRandomNumbers_2D_ArrayOfArrays(double **arr, double row, double col)
{
    /*
    * Fulfilling the array with random numbers 
    * */
    size_t i, j;
    for (i = 0; i < row; i++)
    {
        for (j = 0; j < col; j++)
        {
            arr[i][j] = get_random();
        }
    }
    return arr;
}

void print_2D_ArrayOfArrays(double **arr, int row, int col)
{
    size_t i, j;
    for (i = 0; i < row; i++)
    {
        for (j = 0; j < col; j++)
        {
            printf("%3.3f ", arr[i][j]);
        }
        printf("\n");
    }
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
            x[i] += a[index] * b[j];
        }
    }
    return 0;
}

size_t SequentialMatrixMultiply_ArrayOfArrays(size_t n, double **a, double **b, double **x)
{
    size_t i, j, k;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            x[i][j] = 0;
            for (k = 0; k < n; k++)
            {
                x[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    return 0;
}

void BlockMatrixMultiply_2(size_t N, double **A, double **B, double **C, size_t block)
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            C[i][j] = 0;
        }
    }
    for (int l2 = 0; l2 < N; l2 += block)
    {
        for (int j2 = 0; j2 < N; j2 += block)
        {
#pragma omp parallel for
            for (int i = 0; i < N; i++)
            {
                numberOfThreads = omp_get_num_threads();
                for (int l = l2; l < min(N, l2 + block); l++)
                {
                    for (int j = j2; j < min(N, j2 + block); j++)
                    {
                        C[i][j] += A[i][l] * B[l][j];
                    }
                }
            }
        }
    }
}

void BlockMatrixMultiplication(int n, double **a, double **b, double **c, size_t nb)
{

    // size_t i, j, k, i1, j1, k1;
    // for (i = 0; i < n; i += B)
    //     for (j = 0; j < n; j += B)
    //         for (k = 0; k < n; k += B)
    //             for (i1 = i; i1 < i + B; i1++)
    //                 for (j1 = j; j1 < j + B; j1++)
    //                 {
    //                     double sum = 0;
    //                     for (k1 = k; k1 < k + B; k1++)
    //                     {
    //                         size_t i1_k1 = i1 * n + k1;
    //                         size_t k1_j1 = k1 * n + j1;
    //                         sum += a[i1_k1] * b[k1_j1];
    //                     }
    //                     size_t i1_j1 = i1 * n + j1;
    //                     c[i1_j1] += sum;
    //                 }
    size_t MATRIX_SIZE = n;
    size_t block_size = nb;

    for (int i = 0; i < MATRIX_SIZE; i += block_size)
    {
        for (int j = 0; j < MATRIX_SIZE; j += block_size)
        {
#pragma omp parallel for collapse(2)
            for (int x = 0; x < block_size; ++x)
            {
                for (int y = 0; y < block_size; ++y)
                {
                    for (int k = 0; k < MATRIX_SIZE; ++k)
                    {
                        numberOfThreads = omp_get_num_threads();
#pragma omp critical
                        c[i][j] += a[i][k] * b[k][j];
                    }
                }
            }
        }
    }
}

int main(int argc, char *argv[])
{
    // Global declerations and initilization works
    size_t i, N, NB;
    int world_size = 0;
    double time_start_openmp, time_end_openmp, openmp_exec_time;
    // srand(time(NULL));

    // Parsing program arguments
    if (argc != 5)
    {
        printf("Usage: %s -n <N> -nb <NB>\n", argv[0]);
        return 0;
    }

    for (size_t j = 0; j < 5; j++)
    {
        if (strcmp(argv[j], "-n") == 0)
        {
            N = atoi(argv[j + 1]);
        }
        if (strcmp(argv[j], "-nb") == 0)
        {
            NB = atoi(argv[j + 1]);
        }
    }

    size_t n = N;
    printf("Runing with N=%d and NB=%d\n", n, NB);

    //Array allocations
    double **a = allocarray_2D_ArrayOfArrays(n, n);
    double **b = allocarray_2D_ArrayOfArrays(n, n);
    double **x = allocarray_2D_ArrayOfArrays(n, n);
    double **x_seq = allocarray_2D_ArrayOfArrays(n, n);

    if (a == NULL || b == NULL || x == NULL)
    {
        printf("Allocation failed\n");
        return 0;
    }

    // Array filling
    fullfillArrayWithRandomNumbers_2D_ArrayOfArrays(a, n, n);
    fullfillArrayWithRandomNumbers_2D_ArrayOfArrays(b, n, n);

#ifdef DEBUG
    printf("A:\n");
    print_2D_ArrayOfArrays(a, n, n);
    printf("b:\n");
    print_2D_ArrayOfArrays(b, n, n);
#endif

    SequentialMatrixMultiply_ArrayOfArrays(n, a, b, x_seq);

    // Muliplying matrices
    time_start_openmp = omp_get_wtime();
    BlockMatrixMultiply_2(n, a, b, x, NB);
    time_end_openmp = omp_get_wtime();

#ifdef DEBUG
    printf("\nX:\n");
    print_2D_ArrayOfArrays(x, n, n);
    printf("\nX_seq:\n");
    print_2D_ArrayOfArrays(x_seq, n, n);
#endif

    // Execution time calculations
    openmp_exec_time = time_end_openmp - time_start_openmp;

    // print matrix size, number of processors, number of threads, time, time_openmp, L2 norm of difference of x and xseq (use %.12e while printing norm)
    printf("OPENMP: %d %ld %f\n", n, numberOfThreads, openmp_exec_time, openmp_exec_time);

    // Memory usage calculation (in terms of Gigabytes)
    totalMemUsage = totalMemUsage / (1024 * 1024 * 1024);
    printf("TOTALMEMUSAGE: %zu\n", totalMemUsage);
    return 0;
}
