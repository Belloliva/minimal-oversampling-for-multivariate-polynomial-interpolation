//
//  sampling.h
//  scalar_functions
//
//  Created by abdellah chkifa on 12/24/23.
//

#ifndef mpa_sampling_h
#define mpa_sampling_h

#include <stdio.h>
#include "mpa_monotone_graph.h"
#include "Least_squares_solver.h"


extern size_t num_rejection;


double* iid_uniform_sample_1d (size_t NumSamples);
double* uniform_sample_multi_d (size_t dim);
double** iid_uniform_sample_multi_d (size_t dim, size_t NumSamples);

double* iid_arcsine_sample_1d (size_t NumSamples);
double* arcsine_sample_multi_d (size_t dim);
double** iid_arcsine_sample_multi_d (size_t dim, size_t NumSamples);

double arcsine_pdf_1d (double x);
double arcsine_pdf_multi_d (size_t dim, double* y);

double christoffel_Legendre_kernel_1d(size_t n, double x);
double christoffel_Legendre_pdf_1d (size_t n, double x);
double* iid_christoffel_Legendre_sample_1d (size_t n, size_t N);

void set_Legendre_vector (mpa_monotone_graph* Lambda, size_t activeSize, double *y, gsl_vector* L);

double christoffel_Legendre_kernel_multi_d      (size_t dim, mpa_monotone_graph* Lambda, size_t activeSize, double* y);
double christoffel_Legendre_pdf_multi_d         (size_t dim, mpa_monotone_graph* Lambda, size_t activeSize, double* y);
double skew_christoffel_Legendre_kernel_multi_d (size_t dim, mpa_monotone_graph* Lambda, size_t activeSize, gsl_matrix* M, double* y);
double skew_christoffel_Legendre_pdf_multi_d    (size_t dim, mpa_monotone_graph* Lambda, size_t activeSize, gsl_matrix* M, double* y);

double** iid_christoffel_Legendre_sample_multi_d      (size_t dim, mpa_monotone_graph* Lambda, size_t activeSize, size_t N);
double** iid_skew_christoffel_Legendre_sample_multi_d (size_t dim, mpa_monotone_graph* Lambda, size_t activeSize,
                                                       double eps, double gamma, double alpha, size_t N);

ls_datas* w_uniform_sample_multi_d  (size_t dim, mpa_monotone_graph* Lambda, size_t activeSize, size_t N);
ls_datas* w_arcsine_sample_multi_d(size_t dim, mpa_monotone_graph* Lambda, size_t activeSize, size_t N);
ls_datas* vanilla_CL_sample_multi_d (size_t dim, mpa_monotone_graph* Lambda, size_t activeSize, size_t N);

ls_datas* vbp_CL_sample_multi_d(size_t dim, mpa_monotone_graph* Lambda, size_t activeSize, double eps, double gamma, size_t N);
ls_datas* fbp_CL_sample_multi_d(size_t dim, mpa_monotone_graph* Lambda, size_t activeSize, double delta, double kappa, size_t N);

#endif /* mpa_sampling_h */
