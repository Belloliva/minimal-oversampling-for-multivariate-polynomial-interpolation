//
//  osl_fe_linalg.c
//  ff_tests
//
//  Created by abdellah chkifa on 12/4/23.
//

#include "mpa_linalg.h"
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>


double dot_product (const gsl_vector* left, const gsl_vector* right)
{
    double ans=0;
    gsl_blas_ddot(left, right, &ans);
    return ans;
}

// The matrix A is symmetric
double Quadratic_form (const gsl_matrix* A, const gsl_vector* left, const gsl_vector* right)
{
    gsl_vector* Ar = gsl_vector_calloc(right->size);
    gsl_blas_dsymv(CblasUpper, 1.0, A, right, 0.0, Ar);

    double ans = dot_product(left, Ar);
    gsl_vector_free(Ar);
    return ans;
}

#pragma mark -
#pragma mark functions used in least squares sampling

double matrix_trace (gsl_matrix* A)
{
    double trace=0.0;
    for (int i=0; i<A->size1; i++)
        trace += gsl_matrix_get(A, i, i);
    return trace;
}

gsl_vector* sym_matrix_eigenvalues (const gsl_matrix* A){
    
    size_t n = A->size1;
    gsl_matrix* AA = gsl_matrix_calloc(n, n);
    gsl_matrix_memcpy(AA, A);
    
    gsl_eigen_symm_workspace * W = gsl_eigen_symm_alloc(n);
    gsl_vector *eval = gsl_vector_calloc(n);
    gsl_eigen_symm(AA, eval, W);
    
    gsl_eigen_symm_free(W);
    gsl_matrix_free(AA);
    
    return eval;
}


doublePair sym_matrix_minmax_eigenvalue (const gsl_matrix* A){
    
    gsl_vector* eigenvalues = sym_matrix_eigenvalues(A);
    double lam_min = gsl_vector_min(eigenvalues);
    double lam_max = gsl_vector_max(eigenvalues);
    gsl_vector_free(eigenvalues);
    doublePair ans = {lam_min,lam_max}  ;
    return ans;
}


double sym_matrix_condition_number (const gsl_matrix* A){
    doublePair minmax = sym_matrix_minmax_eigenvalue (A);
    return minmax.second/minmax.first;
}



void matrix_add_multiple_eye (gsl_matrix* A, double l)
{
    for (int i=0; i<A->size1; i++)
        gsl_matrix_set(A, i, i, gsl_matrix_get(A, i, i) + l);
}


void matrix_add_outer_product (gsl_matrix* A, gsl_vector* w, double s)
{
    double wi,wj;
    for (int i=0; i<A->size1; i++){
        wi = gsl_vector_get(w, i);
        for (int j=0; j<A->size2; j++){
            wj = gsl_vector_get(w, j);
            double aij = gsl_matrix_get(A, i, j);
            gsl_matrix_set(A, i, j, aij + s*wi*wj);
        }
    }
}



void set_matrix_resolvant (const gsl_matrix* A, double l, gsl_matrix* RA)
{
    gsl_matrix* B = gsl_matrix_alloc(A->size1, A->size2);
    gsl_matrix_memcpy(B, A);
    matrix_add_multiple_eye (B, -l);
    
    int s;
    gsl_permutation * p = gsl_permutation_alloc (A->size1);

    
    gsl_linalg_LU_decomp (B, p, &s);
    gsl_linalg_LU_invert (B, p, RA);
    
    gsl_matrix_free (B);
    gsl_permutation_free (p);
}

