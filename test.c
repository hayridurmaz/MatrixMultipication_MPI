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

int get_random(int rank)
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

int **allocarray(int row, int col, int which)
{
    int **array1 = malloc(row * sizeof(int *));
    for (int i = 0; i < row; i++)
    {
        array1[i] = malloc(col * sizeof(int));
        printf("col=%d - Size of array1[%d]=%zu (%zu)\n", col, i, sizeof(array1[i]), ((int *)(col * sizeof(int))));
    }

    int *array2[row];
    for (int i = 0; i < row; i++)
    {
        array2[i] = (int *)malloc(col * sizeof(int));
        printf("col=%d - Size of array2[%d]=%zu (%zu)\n", col, i, sizeof(array2[i]), ((int *)(col * sizeof(int))));
    }

    int(*array3)[row] = malloc(sizeof(int[row][col]));

    for (size_t i = 0; i < row; ++i)
        for (size_t j = 0; j < col; ++j)
        {
            array3[i][j] = 1;
            printf("col=%d - Size of array2[%d]=%zu (%zu)\n", col, i, sizeof(array3[i]), ((int *)(col * sizeof(int))));
        }

    int *data = malloc(row * col * sizeof(int));
    int **array4 = malloc(row * sizeof(int *));
    for (int i = 0; i < row; i++)
        array4[i] = &(data[i * col]);

    printf(" in allocarray(array1); %zu - %p mult(%d * %d * %zu)\n", sizeof(array1), array1, row, col, sizeof(int));
    printf(" in allocarray(array2); %zu - %p mult(%d * %d * %zu)\n", sizeof(array2), array2, row, col, sizeof(int));
    printf(" in allocarray(array3); %zu - %p mult(%d * %d * %zu)\n", sizeof(array3), array3, row, col, sizeof(int));
    printf(" in allocarray(array3); %zu - %p mult(%d * %d * %zu)\n", sizeof(array4), array4, row, col, sizeof(int));
    return which ? array1 : array4;
}

int main(int argc, char const *argv[])
{
    /* code */

    int i, j;
    int N = MATSIZE;
    int K, nOverK;
    // Initialize the MPI environment

    // Get the number of processes

    // Get the rank of the process
    int taskid;

    // Get the name of the processor
    int name_len;

    K = 4;
    nOverK = N / K;

    int **arr = allocarray(nOverK, N, 1);
    printf("after allocarray:: %d\n", sizeof(arr));
    srand(time(NULL) + taskid);
    /*
        * Fulfilling the array with random numbers 
        * */
    for (i = 0; i < nOverK; i++)
    {
        for (j = 0; j < N; j++)
        {
            arr[i][j] = get_random(1);
        }
    }
    print_arr(arr, nOverK, N);
    free(arr);

    arr = allocarray(nOverK, N, 0);
    printf("after allocarray:: %d\n", sizeof(arr));
    srand(time(NULL) + taskid);
    /*
        * Fulfilling the array with random numbers 
        * */
    for (i = 0; i < nOverK; i++)
    {
        for (j = 0; j < N; j++)
        {
            arr[i][j] = get_random(1);
        }
    }
    print_arr(arr, nOverK, N);
    free(arr);

    return 0;
}
