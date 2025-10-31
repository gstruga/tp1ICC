#include <stdio.h>
#include <stdlib.h> /* for exit e random/srandom */
#include <string.h>
#include <float.h>
#include <math.h>

#include "utils.h"
#include "sislin.h"

#define LIMITEI(k) (((int)(ceil(k / 2.0))) - (1))
#define LIMITEJ(k) ((int)(floor(k / 2.0)))

#define SEM_PRE_CONDICIONADOR -1.0
#define JACOBI 0.0
#define GAUSS_SEIDEL 1.0

/* ----- FUNCOES AUXILIARES ----- */

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

/*faz o vetor trasposto vezes o vetor*/
static inline real_t produtoInternoDeVetores(real_t *vetor1, real_t *vetor2, int tam)
{
  real_t acomulador = 0.0;
  for (int i = 0; i < tam; i++)
  {
    acomulador += vetor1[i] * vetor2[i];
  }
  return acomulador;
}

/*faz o vetor transposto vezes uma matriz vezes o vetor*/
real_t vetorTxMatrizxVetor(real_t *vetor, real_t **Matriz, int tam, int k)
{
  int limitei = LIMITEI(k);
  int limitej = LIMITEJ(k);
  real_t acumuladores[tam], acumulador = 0;
  memset(acumuladores, 0, tam * sizeof(real_t));

  for (int i = 0; i < tam; i++)
  {
    int j = 0;
    while (j < tam)
    {
      if (i > limitei + j)
      {
        j++;
      }
      else if (i < j - limitej)
      {
        break;
      }
      else
      {
        acumuladores[i] += vetor[j] * Matriz[i][j];
        j++;
      }
    }
  }
  for (int i = 0; i < tam; i++)
  {
    acumulador += vetor[i] * acumuladores[i];
  }
  return acumulador;
}

/* faz matriz vezes vetor */
void matrizXvetor(real_t **matriz, real_t *vetor, real_t *vetorDestino, int tam, int k)
{
  int limitei = LIMITEI(k);
  int limitej = LIMITEJ(k);

  memset(vetorDestino, 0, tam * sizeof(real_t));

  for (int i = 0; i < tam; i++)
  {
    int j = 0;
    while (j < tam)
    {
      if (i > limitei + j)
      {
        j++;
      }
      else if (i < j - limitej)
      {
        break;
      }
      else
      {
        vetorDestino[i] += vetor[j] * matriz[i][j];
        j++;
      }
    }
  }
}

/*calcula a norma maxima do erro*/
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

/* soma a matriz com sua transposta e guarda em resultado */
void somaComTransposta(real_t **resultado, real_t **matriz, int n)
{
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
  int limitei = LIMITEI(k);
  int limitej = LIMITEJ(k);

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
  real_t **AT = aloca_matriz(n);
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
    {
      AT[i][j] = A[j][i];
    }
  }

  multiplicaPeloPreCondicionador(ASP, A, AT, bsp, b, n, NULL);
  /*int limitej = LIMITEJ(k);
  for (int i = 0; i < n; i++)
  {
    ASP[i][i] *= 0.5;
    for (int j = i + 1; j <= i + limitej; j++)
    {
      if (j >= n)
        break;
      ASP[i][j] *= 0.5;
    }

    for (int j = i - 1; j <= i - limitej; j--)
    {
      if (j < 0)
        break;
      ASP[i][j] *= 0.5;
    }
  }

  memcpy(bsp, b, n * sizeof(real_t));*/

  for (int i = 0; i < n; i++)
  {
    free(AT[i]);
  }
  free(AT);

  *tempo = timestamp() - *tempo;
}

void geraDLU(real_t **A, int n, int k,
             real_t *D, real_t **L, real_t **U, rtime_t *tempo)
{
  *tempo = timestamp();
  int limitej = LIMITEJ(k);
  for (int i = 0; i < n; i++)
  {
    D[i] = A[i][i];
    for (int j = i + 1; j <= i + limitej; j++)
    {
      if (j >= n)
        break;

      U[i][j] = A[i][j];
    }
    for (int j = i - 1; j >= i - limitej; j--)
    {
      if (j < 0)
        break;

      L[i][j] = A[i][j];
    }
  }

  *tempo = timestamp() - *tempo;
}

/**
 * Devolve matriz M⁻¹
 *
 */
void geraPreCond(real_t *D, real_t **L, real_t **U, real_t w, int n, int k,
                 real_t **M, rtime_t *tempo)
{
  *tempo = timestamp();

  if (w == SEM_PRE_CONDICIONADOR)
  {
    for (int i = 0; i < n; i++)
    {
      M[i][i] = 1.0;
    }
    return;
  }
  else if (w == JACOBI)
  {
    for (int i = 0; i < n; i++)
    {
      M[i][i] = 1.0 / D[i];
    }
    return;
  }
  else
  {
    perror("OMEGA INVALIDO\n");
    exit(-1);
  }

  *tempo = timestamp() - *tempo;
}

/*preenche o vetor residuo e retorna norma L2 do residuo*/
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

/*aplica o gradiente conjugado e retorna a norma da solucao com a iteracao anterior*/
real_t gradienteConjugado(real_t **A, real_t *b, real_t *X, real_t **M,
                          int n, int k, rtime_t *tempo, int max_it, real_t epsilon, rtime_t *tempoResiduo, int *totalIt)
{
  real_t *residuo, *p, *Apk, *XAntigo, *zk, betaK;
  real_t norma;
  residuo = malloc(sizeof(real_t) * n);
  if (!residuo)
  {
    perror("ERRO AO ALOCAR VETOR\n");
    exit(-1);
  }

  XAntigo = calloc(sizeof(real_t), n);
  if (!XAntigo)
  {
    perror("ERRO AO ALOCAR VETOR\n");
    exit(-1);
  }

  zk = calloc(sizeof(real_t), n);
  if (!zk)
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

  Apk = malloc(sizeof(real_t) * n);
  if (!Apk)
  {
    perror("ERRO AO ALOCAR VETOR\n");
    exit(-1);
  }

  *tempo = timestamp();

  *totalIt = 0;
  if (calcResiduoSL(A, b, X, n, k, residuo, tempo) - epsilon <= DBL_EPSILON)
    return 0.0;

  matrizXvetor(M, residuo, zk, n, k);
  memcpy(p, zk, n * sizeof(real_t));
  for (int i = 0; i < max_it; i++)
  {
    real_t numerador = produtoInternoDeVetores(residuo, zk, n);
    matrizXvetor(A, p, Apk, n, k);
    real_t ak = numerador / vetorTxMatrizxVetor(p, A, n, k);

    /* matrizXvetor(A, p, Apk, n, k); */

    memcpy(XAntigo, X, n * sizeof(real_t));
    for (int j = 0; j < n; j++)
    {
      X[j] += ak * p[j];
      residuo[j] -= ak * Apk[j];
    }

    (*totalIt)++;

    *tempoResiduo = timestamp();
    real_t residuoNovo = produtoInternoDeVetores(residuo, residuo, n);
    *tempoResiduo = timestamp() - *tempoResiduo;
    norma = norma_erro_maximo(n, XAntigo, X);
    if (norma - epsilon < DBL_EPSILON)
    {
      break;
    }
    matrizXvetor(M, residuo, zk, n, k);

    betaK = produtoInternoDeVetores(residuo, zk, n) / numerador;

    for (int j = 0; j < n; j++)
    {
      p[j] = zk[j] + betaK * p[j];
    }
  }

  *tempo = timestamp() - *tempo;

  free(XAntigo);

  free(residuo);

  free(Apk);

  free(p);

  free(zk);

  return norma;
}

/* multiplica a matriz de coeficientes e os termos independentes pelo precondicionador
 * e guarda em matrizNova e em vetorNovo */
void multiplicaPeloPreCondicionador(real_t **matrizNova, real_t **A, real_t **M, real_t *vetorNovo, real_t *B, int n, rtime_t *tempo)
{
  //*tempo = timestamp();
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
    {
      matrizNova[i][j] = 0.0;
      for (int k = 0; k < n; k++)
      {
        matrizNova[i][j] += M[i][k] * A[k][j];
      }
    }
  }

  for (int i = 0; i < n; i++)
  {
    vetorNovo[i] = 0.0;
    for (int k = 0; k < n; k++)
    {
      vetorNovo[i] += M[i][k] * B[k];
    }
  }

  //*tempo = timestamp() - *tempo;
}