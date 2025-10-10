#include <stdio.h>
#include <stdlib.h> /* for exit e random/srandom */
#include <string.h>
#include <float.h>
#include <math.h>

#include "utils.h"
#include "sislin.h"

/* ----- FUNCOES AUXILIARES ----- */

real_t norma_erro_maximo(int n, double *x_k, double *x_k_1)
{
  real_t max_diferenca = 0.0;

  for (int i = 0; i < n; i++)
  {
    real_t diff = fabs(x_k[i] - x_k_1[i]);
    if (diff > max_diferenca)
    {
      max_diferenca = diff;
    }
  }

  return max_diferenca;
}

void somaComTransposta(real_t** resultado, real_t** matriz, int n){
  if (!resultado || !matriz)
  {
    perror("UM DOS ARGUMENTOS DA FUNCAO EH NULL\n");
    exit(-1);
  }
  
  
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
    {
      resultado[i][j] = matriz[i][j] + matriz[j][i];
    }
    
  }
  
}

/* ----- FUNCOES AUXILIARES ----- */

static inline real_t generateRandomA(unsigned int i, unsigned int j, unsigned int k);
static inline real_t generateRandomB(unsigned int k);
static real_t invRandMax = 1.0 / (real_t)RAND_MAX;
/**
 * Função que gera os coeficientes de um sistema linear k-diagonal
 * @param i,j coordenadas do elemento a ser calculado (0<=i,j<n)
 * @param k numero de diagonais da matriz A
 */
static inline real_t generateRandomA(unsigned int i, unsigned int j, unsigned int k)
{
  return ((i == j) ? (real_t)(k << 1) : 1.0) * (real_t)random() * invRandMax;
}

/**
 * Função que gera os termos independentes de um sistema linear k-diagonal
 * @param k numero de diagonais da matriz A
 */
static inline real_t generateRandomB(unsigned int k)
{
  return (real_t)(k << 2) * (real_t)random() * invRandMax;
}

/* Cria matriz 'A' k-diagonal e Termos independentes B */
void criaKDiagonal(int n, int k, real_t **A, real_t *B)
{
  int limitei = (int)ceil(k / 2.0) - 1;
  int limitej = (int)floor(k / 2.0);

  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
    {
      if (i > limitei + j || i < j - limitej)
      {
        A[i][j] = 0.0;
      }
      else
      {
        A[i][j] = generateRandomA(i, j, k);
      }
    }
  }

  for (int i = 0; i < n; i++)
    B[i] = generateRandomB(k);
}

/* Gera matriz simetrica positiva */
void genSimetricaPositiva(real_t **A, real_t *b, int n, int k,
                          real_t **ASP, real_t *bsp, rtime_t *tempo)
{
  *tempo = timestamp();
  somaComTransposta(ASP, A, n);

  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
    {
      ASP[i][j] *= 0.5;
    }
    
  }

   for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            printf("%f ", ASP[i][j]);
        }

        printf("\n");
        
    }
  
  *tempo = timestamp() - *tempo;
}

void geraDLU(real_t *A, int n, int k,
             real_t **D, real_t **L, real_t **U, rtime_t *tempo)
{
  *tempo = timestamp();

  *tempo = timestamp() - *tempo;
}

/**
 * Devolve matriz M⁻¹
 *
 */
void geraPreCond(real_t *D, real_t *L, real_t *U, real_t w, int n, int k,
                 real_t **M, rtime_t *tempo)
{
  *tempo = timestamp();

  *tempo = timestamp() - *tempo;
}

real_t calcResiduoSL(real_t **A, real_t *b, real_t *X,
                     int n, int k, real_t *residuo, rtime_t *tempo)
{
  *tempo = timestamp();

  real_t soma = 0.0;

  for (int i = 0; i < n; i++)
  {
    real_t Ax_i = 0.0;

    for (int j = 0; j < n; j++)
    {
      Ax_i += A[i][j] * X[j];
    }

    real_t ri = b[i] - Ax_i;
    residuo[i] = ri;
    soma += ri * ri;
  }

  *tempo = timestamp() - *tempo;
  return sqrt(soma);
}

void resolveSemPreCondicionador(real_t **A, real_t *b, real_t *X,
                                int n, int k, rtime_t *tempo, int max_it, real_t epsilon)
{
  real_t **ASP, *BSP, *residuo, *p;
  ASP = (real_t **)malloc(sizeof(real_t *) * n);

  if (!ASP)
  {
    perror("ERRO AO ALOCAR MATRIZ\n");
    exit(-1);
  }

  BSP = malloc(sizeof(real_t) * n);
  if (!BSP)
  {
    perror("ERRO AO ALOCAR VETOR\n");
    exit(-1);
  }

  for (int i = 0; i < n; i++)
  {
    ASP[i] = malloc(sizeof(real_t) * n);
    if (!ASP[i])
    {
      perror("ERRO AO ALOCAR VETOR\n");
      exit(-1);
    }
  }

  residuo = malloc(sizeof(real_t) * n);
  if (!residuo)
  {
    perror("ERRO AO ALOCAR VETOR\n");
    exit(-1);
  }

  p = malloc(sizeof(real_t) * n);
  if (!p)
  {
    perror("ERRO AO ALOCAR VETOR\n");
    exit(-1);
  }

  *tempo = timestamp();

  genSimetricaPositiva(A, b, n, k, ASP, BSP, tempo);
  if(calcResiduoSL(ASP, b, X, n, k, residuo, tempo) - epsilon <= DBL_EPSILON)
    return;

  memcpy(p, residuo, n * sizeof(real_t));
  for (int i = 0; i < max_it; i++)
  {
    real_t ak = multiplicaComTransposta(residuo, n) / vetorMultComMatriz(p, A, n);

    for (int j = 0; j < n; j++)
    {
      X[j] += ak * p[j];
    }
    
  }
  
  *tempo = timestamp() - *tempo;

  free(residuo);
  free(BSP);
 
  for (int i = 0; i < n; i++)
  {
    free(ASP);
  }

  free(ASP);
  free(p);
}