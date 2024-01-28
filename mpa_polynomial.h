//
//  mpa_polynomial.h
//  thesis_tests
//
//  Created by abdellah chkifa on 11/3/23.
//  Copyright © 2023 chkifa abdellah. All rights reserved.
//

#ifndef mpa_polynomial_h
#define mpa_polynomial_h

#include <stdio.h>
#include "mpa_multi_index.h"


double univariate_monomial              (size_t k, double x, void* unused);
double univariate_Chebychev_polynomial  (size_t k, double x, void* unused);
double univariate_Legendre_polynomial_P (size_t k, double x, void* unused);
double univariate_Legendre_polynomial_L (size_t k, double x, void* unused);

double multivariate_monomial              (const mpa_multi_index* nu, const double* xvalue, void* unused);
double multivariate_Chebychev_polynomial  (const mpa_multi_index* nu, const double* xvalue, void* unused);
double multivariate_Legendre_polynomial_P (const mpa_multi_index* nu, const double* xvalue, void* unused);
double multivariate_Legendre_polynomial_L (const mpa_multi_index* nu, const double* xvalue, void* unused);

#endif /* mpa_polynomial_h */
