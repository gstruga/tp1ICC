#include <stdio.h>
#include "sislin.h"
#include "utils.h"

/*libera toda memoria da heap*/
void free_tudo(real_t **A, real_t **ASP, real_t **L, real_t **U, real_t **M, real_t *BSP,
               real_t *B, real_t *X, real_t *residuo, real_t *D, int n)
{
    for (int i = 0; i < n; i++)
    {
        free(A[i]);
        free(ASP[i]);
        free(L[i]);
        free(U[i]);
        free(M[i]);
    }
    free(A);
    free(ASP);
    free(BSP);
    free(B);
    free(X);
    free(residuo);
    free(D);
    free(L);
    free(U);
    free(M);
}

/*aloca uma matriz na heap*/
static inline real_t **aloca_matriz(int tam)
{
    real_t **aux = (real_t **)calloc(sizeof(real_t *), tam);

    if (!aux)
    {
        perror("ERRO AO ALOCAR MATRIZ\n");
        exit(-1);
    }
    for (int i = 0; i < tam; i++)
    {
        aux[i] = (real_t *)calloc(sizeof(real_t), tam);
        if (!aux[i])
        {
            perror("ERRO AO ALOCAR VETOR\n");
            exit(-1);
        }
    }

    return (aux);
}

/*aloca um vetor na heap*/
static inline real_t *aloca_vetor(int tam)
{
    real_t *aux = (real_t *)calloc(sizeof(real_t), tam);

    if (!aux)
    {
        perror("ERRO AO ALOCAR VETOR\n");
        exit(-1);
    }

    return aux;
}

void prnMat(real_t **mat, real_t *vet, int tam)
{
    for (int i = 0; i < tam; i++)
    {
        for (int j = 0; j < tam; j++)
        {
            printf("%lf ", mat[i][j]);
        }
        printf("%lf\n", vet[i]);
    }
}

int main()
{
    srandom(20252);
    int n, k, maxit, totalIt;
    real_t w = 0.0;
    real_t **A, *B, **ASP, *BSP, *residuo;
    real_t *D, **L, **U, **M, e;
    rtime_t tempoPC, tempoIter, tempoResiduo;

    scanf("%d %d %lf %d %lf", &n, &k, &w, &maxit, &e);
    if (w != -1.0 && w != 0.0)
    {
        printf("NAO IMPLEMENTEI PRECONDICIONADOR DE GAUSS-SEIDEL E SSOR\n");
        return -1;
    }

    A = aloca_matriz(n);
    B = aloca_vetor(n);
    residuo = aloca_vetor(n);
    ASP = aloca_matriz(n);
    BSP = aloca_vetor(n);

    D = aloca_vetor(n);
    L = aloca_matriz(n);
    U = aloca_matriz(n);
    M = aloca_matriz(n);

    real_t *X = aloca_vetor(n);

    criaKDiagonal(n, k, A, B);

    // prnMat(A, B, n);

    double tempo;

    tempoPC = timestamp();
    genSimetricaPositiva(A, B, n, k, ASP, BSP, &tempo);
    // printf("\n\n");
    // prnMat(ASP, BSP, n);
    geraDLU(ASP, n, k, D, L, U, &tempo);
    geraPreCond(D, L, U, w, n, k, M, &tempo);
    tempoPC = timestamp() - tempoPC;

    tempoIter = timestamp();
    real_t norma = gradienteConjugado(ASP, BSP, X, M, n, k * 2 - 1, &tempo, maxit, e, &tempoResiduo, &totalIt);
    tempoIter = timestamp() - tempoIter;

    tempoIter = tempoIter / totalIt;
    printf("%d\n", n);
    for (int i = 0; i < n; i++)
    {
        printf("%.16g ", X[i]);
    }
    printf("\n");

    printf("%.8g\n", norma);

    real_t res = calcResiduoSL(A, B, X, n, k, residuo, &tempo);

    printf("%.8g\n", res);

    printf("%.8g\n", tempoPC);

    printf("%.8g\n", tempoIter);

    printf("%.8g\n", tempoResiduo);

    free_tudo(A, ASP, L, U, M, BSP, B, X, residuo, D, n);

    return 0;
}