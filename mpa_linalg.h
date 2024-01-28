//
//  mpa_linalg.h
//  
//
//  Created by abdellah chkifa on 12/4/23.
//

#ifndef mpa_linalg_h
#define mpa_linalg_h

#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>


typedef struct{
    double first;
    double second;
}doublePair;

double dot_product    (const gsl_vector* left, const gsl_vector* right);
double Quadratic_form (const gsl_matrix* A, const gsl_vector* left, const gsl_vector* right);


double matrix_trace (gsl_matrix* A);
doublePair sym_matrix_minmax_eigenvalue (const gsl_matrix* A);
double sym_matrix_condition_number      (const gsl_matrix* A);
gsl_vector* sym_matrix_eigenvalues      (const gsl_matrix* A);

void matrix_add_multiple_eye (gsl_matrix* A, double l);
void matrix_add_outer_product (gsl_matrix* A, gsl_vector* w, double s);

void set_matrix_resolvant (const gsl_matrix* A, double l, gsl_matrix* RA);


#endif /* mpa_linalg_h */
