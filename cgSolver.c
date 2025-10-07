#include <stdio.h>
#include "sislin.h"
#include "utils.h"

int main()
{
    srand(20252);
    int n = 10;
    real_t **A, *B;

    A = (real_t **) malloc(sizeof(real_t *) * n);
    B = malloc(sizeof(real_t) * n);

    for (int i = 0; i < n; i++)
    {
        A[i] = malloc(sizeof(real_t) * n);
    }
    
    criaKDiagonal(n, 3, A, B);

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            printf("%lld ", A[i][j]);
        }

        printf("%lld\n", B[i]);
        
    }
    
    return 0;
}