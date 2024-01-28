//
//  ls_datas.h
//  scalar_functions
//
//  Created by abdellah chkifa on 12/23/23.
//

#ifndef Least_squares_solver_h
#define Least_squares_solver_h

#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_matrix.h>
#include "mpa_monotone_graph.h"
#include "mpa_polynomial.h"

typedef struct ls_datas ls_datas;

struct ls_datas{
    size_t dim; // the underlying dimension of inputs to the target function
    size_t numPoints; // the number of points used in the least squares
    double** points; // the points used in in the least squares
    double* weights;
    gsl_matrix* designMatrix;
    gsl_matrix* gramMatrix;
};

ls_datas* ls_datas_alloc(size_t dim, size_t numPoints);
void ls_datas_free(ls_datas* datas);

void ls_datas_set_points (ls_datas* datas, double** points);
void ls_datas_set_weights(ls_datas* datas, double* weights);

void ls_datas_set_designMatrix (ls_datas* datas, gsl_matrix* D);
void ls_datas_set_gramMatrix   (ls_datas* datas, gsl_matrix* G);

void ls_datas_compute_designMatrix (ls_datas* datas, mpa_monotone_graph* Lambda, size_t active_size);
void ls_datas_compute_gramMatrix   (ls_datas* datas);

void set_ls_solution (ls_datas* datas, mpa_monotone_graph* Lambda, size_t active_size, query_solver qsolver);

#endif /* Least_squares_solver_h */
