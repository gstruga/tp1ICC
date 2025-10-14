#include <stdio.h>
#include "sislin.h"
#include "utils.h"

int main()
{
    srandom(20252);
    int n = 8;
    int k = 5;
    real_t **A, *B, *residuo;

    A = (real_t **)malloc(sizeof(real_t *) * n);
    B = malloc(sizeof(real_t) * n);
    residuo = malloc(sizeof(real_t) * n);

    for (int i = 0; i < n; i++)
    {
        A[i] = malloc(sizeof(real_t) * n);
    }

    criaKDiagonal(n, k, A, B);

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            printf("%f ", A[i][j]);
        }

        printf("%f\n", B[i]);
    }

    double tempo;
    real_t *X = calloc(n, sizeof(real_t));

    resolveSemPreCondicionador(A, B, X, n, k, &tempo, 20, __DBL_EPSILON__);

    for (int i = 0; i < n; i++)
    {
        printf("%f\n", X[i]);
    }

    real_t res = calcResiduoSL(A, B, X, n, k, residuo, &tempo);

    printf("%f\n", res);

    for (int i = 0; i < n; i++)
    {
        free(A[i]);
    }
    free(A);
    free(B);
    free(X);
    free(residuo);

    return 0;
}