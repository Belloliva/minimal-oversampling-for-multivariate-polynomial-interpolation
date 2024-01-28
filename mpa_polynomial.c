//
//  mpa_polynomial.c
//  thesis_tests
//
//  Created by abdellah chkifa on 11/3/23.
//  Copyright © 2023 chkifa abdellah. All rights reserved.
//

#include <math.h>
#include "mpa_polynomial.h"


/*
 Evaluate the polynomial W_k = x^k
 @param k : the index of the monomial
 @param x : the point where the polynomial is evaluated
 @return the value x^k

*/
double univariate_monomial(size_t k, double x, void* unused)
{
    if (k==0) return 1.;
    else return pow(x,k);
}

/*
 Evaluate the tensorized monomial  x_1^{nu_1} x_d^{nu_d}
 @param nu : the multi-index of the  multi-variate monomial
 @param x  : the vector where the polynomial is evaluated, of the same length as the multi-index
 @return the value
 */
double multivariate_monomial(const mpa_multi_index* nu, const double* xvalue, void* unused){
    double prod=1;
    for (size_t j=0; j<nu->size; j++){
        size_t i    = nu->data[j].index;
        size_t nu_i = nu->data[j].coordinate;
        prod*= univariate_monomial(nu_i, xvalue[i], unused);
    }
    return prod;
}

/*
 Evaluate the univariate Chebychev polynomial T_k(x)
 @param k : the index of the polynomial
 @param x : the point where the polynomial is evaluated
 @return the value
*/

double univariate_Chebychev_polynomial(size_t k, double x, void* unused)
{
    if (k==0) return 1.;
    else return cos(k*acos(x));
}

/*
 Evaluate the multivariate Chebychev polynomial  T_nu (x) = T_{nu_1}(x_1)...T_{nu_d}(x_d)
 @param nu : the multi-index of the multi-variate polynomial
 @param x  : the vector where the polynomial is evaluated, of the same length as the multi-index
 @return the value
 */
double multivariate_Chebychev_polynomial(const mpa_multi_index* nu, const double* xvalue, void* unused)
{
    double prod=1;
    for (size_t j=0; j<nu->size; j++){
        size_t i    = nu->data[j].index;
        size_t nu_i = nu->data[j].coordinate;
        prod*= univariate_Chebychev_polynomial(nu_i, xvalue[i], unused);
    }
    return prod;
}


/*
 Evaluate the univariate Legendre polynomial P_k(x) normalized with Pk(1)=1
 @param k : the index of the polynomial
 @param x : the point where the polynomial is evaluated
 @return the value
*/

double univariate_Legendre_polynomial_P(size_t k, double x, void* unused)
{
    if(k==0) return 1.;
    if(k==1) return x ;
    double P0 = 1, P1=x, Px=0.0;
    for (int j=2; j<=k; j++) {
        Px = ((2*j-1)*x*P1 - (j-1)*P0)/j;
        P0 = P1;
        P1 = Px;
    }
    return Px;
}

/*
 Evaluate the multivariate Legendre polynomial  P_nu (x) = P_{nu_1}(x_1)...P_{nu_d}(x_d)
 @param nu : the multi-index of the multi-variate polynomial
 @param x  : the vector where the polynomial is evaluated, of the same length as the multi-index
 @return the value
 */
double multivariate_Legendre_polynomial_P(const mpa_multi_index* nu, const double* xvalue, void* unused)
{
    double prod=1;
    for (size_t j=0; j<nu->size; j++){
        size_t i    = nu->data[j].index;
        size_t nu_i = nu->data[j].coordinate;
        prod*= univariate_Legendre_polynomial_P(nu_i, xvalue[i], unused);
    }
    return prod;
}

/*
 Evaluate the univariate Legendre polynomial L_k(x) normalized with L2 norm equal 1,
 w.r.t. measure dx/2 over [-1,1]
 @param k : the index of the polynomial
 @param x : the point where the polynomial is evaluated
 @return the value
*/

double univariate_Legendre_polynomial_L(size_t k, double x, void* unused)
{
    return sqrt(2*k+1)*univariate_Legendre_polynomial_P(k,x, unused);
}

/*
 Evaluate the multivariate Legendre polynomial  L_nu (x) = L_{nu_1}(x_1)...L_{nu_d}(x_d)
 @param nu : the multi-index of the multi-variate polynomial
 @param x  : the vector where the polynomial is evaluated, of the same length as the multi-index
 @return the value
 */
double multivariate_Legendre_polynomial_L(const mpa_multi_index* nu, const double* xvalue, void* unused){
    double prod=1;
    for (size_t j=0; j<nu->size; j++){
        size_t i    = nu->data[j].index;
        size_t nu_i = nu->data[j].coordinate;
        prod*= univariate_Legendre_polynomial_L(nu_i, xvalue[i], unused);
    }
    return prod;
}
