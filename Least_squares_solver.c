//
//  Least_squares_solver.c
//  scalar_functions
//
//  Created by abdellah chkifa on 12/23/23.
//

#include "Least_squares_solver.h"
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

ls_datas* ls_datas_alloc(size_t dim, size_t numPoints)
{
    ls_datas* datas = (ls_datas*) malloc(sizeof(ls_datas));
    datas->dim = dim;
    datas->numPoints = numPoints;
    return datas;
}

void ls_datas_free(ls_datas* datas)
{
    for (int i=0; i<datas->numPoints; i++)
        free(datas->points[i]);
    free(datas->points);
    free(datas->weights);

    gsl_matrix_free(datas->designMatrix);
    gsl_matrix_free(datas->gramMatrix);
}

void ls_datas_set_points (ls_datas* datas, double** points){ datas->points = points;}
void ls_datas_set_weights(ls_datas* datas, double* weights){datas->weights = weights;}
void ls_datas_set_designMatrix (ls_datas* datas, gsl_matrix* D){datas->designMatrix= D;}
void ls_datas_set_gramMatrix   (ls_datas* datas, gsl_matrix* G){datas->gramMatrix  = G;}

void ls_datas_compute_designMatrix (ls_datas* datas, mpa_monotone_graph* Lambda, size_t active_size)
{
    datas->designMatrix = gsl_matrix_calloc(datas->numPoints, active_size);

    for (int j=0; j< active_size; j++){
        mpa_multi_index* nu = Lambda->nodes[j]->miPtr;
        for (int i=0; i<datas->numPoints; i++){
            double* y = datas->points[i];
            double Ly = multivariate_Legendre_polynomial_L(nu, y, NULL);
            gsl_matrix_set(datas->designMatrix, i, j, Ly);
        }
    }
}

void ls_datas_compute_gramMatrix (ls_datas* datas){

    gsl_matrix * A = datas->designMatrix;
    size_t n = A->size2;
    datas->gramMatrix = gsl_matrix_calloc(n, n);
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1., A, A, 0., datas->gramMatrix);
}


void set_ls_solution(ls_datas* datas, mpa_monotone_graph* Lambda, size_t active_size,
                      query_solver qsolver)
{
    gsl_matrix * G = datas->gramMatrix;
    double** x = datas->points;
    double* s  = datas->weights;
    size_t numPoints = datas->numPoints;
    
    gsl_vector* vals = gsl_vector_calloc(numPoints);
    for (int i=0; i<numPoints; i++){
        double fy=0;
        double* y = x[i];
        qsolver.evaluate (y, qsolver.dim, qsolver.datas, &fy);
        gsl_vector_set(vals, i, fy*s[i]);
    }

    gsl_vector* rhs = gsl_vector_calloc(active_size);
        
    gsl_matrix* A = datas->designMatrix;
    gsl_blas_dgemv(CblasTrans, 1., A, vals, 0, rhs);
    
    
    gsl_linalg_cholesky_decomp(G);
    gsl_vector * result=gsl_vector_calloc(active_size);
    gsl_linalg_cholesky_solve(G, rhs, result);
    
    for (int i=0; i<active_size; i++) {
        Lambda->nodes[i]->info = (double*) malloc(sizeof(double));
        double* info = Lambda->nodes[i]->info;
        *info = gsl_vector_get(result, i);
        //printf("coef at %d is %f \n",i,gsl_vector_get(result, i));
    }
    gsl_vector_free(rhs);
    gsl_vector_free(result);
}



