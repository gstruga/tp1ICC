#ifndef __SISLIN_H__
#define __SISLIN_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "utils.h"
#include "sislin.h"

void criaKDiagonal(int n, int k, double **A, double *B);

void genSimetricaPositiva(double **A, double *b, int n, int k, double **ASP, double *bsp, double *tempo);
void geraDLU(double **A, int n, int k, double *D, double **L, double **U, double *tempo);
void geraPreCond(double *D, double **L, double **U, double w, int n, int k, double **M, double *tempo);
real_t calcResiduoSL(real_t **A, real_t *b, real_t *X, int n, int k, real_t *residuo, rtime_t *tempo);
real_t gradienteConjugado(real_t **A, real_t *b, real_t *X, real_t **M, int n, int k, rtime_t *tempo, int max_it, real_t epsilon, rtime_t *tempoResiduo, int *totalIt);

#endif // __SISLIN_H__
