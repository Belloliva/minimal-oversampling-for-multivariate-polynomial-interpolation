//
//  map_sampling.c
//  mpa
//
//  Created by abdellah chkifa on 12/24/23.
//

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>

#include "mpa_sampling.h"
#include "mpa_polynomial.h"
#include "mpa_linalg.h"


#pragma mark -
#pragma mark sampling from 1d and multi_d uniform measure over [-1,1]

double uniform_sample_1d (void)
{
    return 2.0*rand()/RAND_MAX -1.0;
}

double* iid_uniform_sample_1d (size_t N)
{
    double* v = (double*) calloc(N, sizeof(double)) ;
    for (int i=0; i<N; i++)
        v[i] = 2.*rand()/RAND_MAX - 1.;
    return v;
}


double* uniform_sample_multi_d (size_t dim)
{
    double* v = (double*) calloc(dim, sizeof(double));
    for (int j=0; j<dim; j++)
        v[j] = 2.*rand()/RAND_MAX - 1.;
    return v;
}


double** iid_uniform_sample_multi_d (size_t dim, size_t N)
{
    double** v = (double**) calloc(N, sizeof(double*));
    for (int i=0; i<N; i++){
        v[i] = (double*) calloc(dim, sizeof(double)) ;
        for (int j=0; j<dim; j++)
            v[i][j] = 2.0*rand()/RAND_MAX-1;
    }
    return v;
}


double uniform_pdf_1d (double x)
{
    return 0.5;
}

#pragma mark -
#pragma mark sampling from 1d and multi_d arcsine measure over [-1,1]

double arcsine_sample_1d (void)
{
    return cos(M_PI*rand()/RAND_MAX);
}

double* iid_arcsine_sample_1d (size_t N)
{
    double* v = (double*) calloc(N, sizeof(double)) ;
    for (int i=0; i<N; i++)
        v[i] = cos(M_PI*rand()/RAND_MAX);
    return v;
}

double* arcsine_sample_multi_d (size_t dim)
{
    double* v = (double*) calloc(dim, sizeof(double));
    for (int j=0; j<dim; j++)
        v[j] = cos(M_PI*rand()/RAND_MAX);
    return v;
}

double** iid_arcsine_sample_multi_d (size_t dim, size_t N)
{
    double** v = (double**) calloc(N, sizeof(double*));
    for (int i=0; i<N; i++){
        v[i] = (double*) calloc(dim, sizeof(double)) ;
        for (int j=0; j<dim; j++)
            v[i][j] = cos(M_PI*rand()/RAND_MAX);
    }
    return v;
}

double arcsine_pdf_1d (double x)
{
    return 1.0/(M_PI* sqrt(1-x*x));
}

double arcsine_pdf_multi_d (size_t dim, double* y)
{
    double prod = 1.;
    for (int i=0; i<dim; i++)
        prod *= arcsine_pdf_1d(y[i]);
    return prod;
}

#pragma mark -
#pragma mark sampling from 1d and multi_d squared Legendre polynomial based measure

double squared_L_pdf_1d (size_t k, double x)
{
    double Lk_x = univariate_Legendre_polynomial_L(k, x, NULL);
    return 0.5* Lk_x*Lk_x;
}


double squared_L_pdf_multi_d (size_t dim,mpa_multi_index* nu, double* x)
{
    double Lnu_x = multivariate_Legendre_polynomial_L(nu, x, NULL);
    return Lnu_x*Lnu_x/pow(2.0, dim);
}

// sample from density |L_k(x)|^2 dx/2 for k>=1
// using rejection method based on Berstein inequality
// for Legendre polynomials
double squared_L_sample_1d (size_t k)
{
    double c = (2.0*k+1.0)/k;
    double x, u, fx, gx;
    int i=0;
    //printf("---------------------\n");
    while (1) {
        u  = 1.0*rand()/RAND_MAX;
        x =  arcsine_sample_1d();
        fx = squared_L_pdf_1d(k, x);
        gx = arcsine_pdf_1d(x);
        if (fx > c *gx*u)
            break;
        else{
            //printf("------%d\n",i);
            i+=1;
        }
    }
    return x;
}

// sample from density |L_nu(x)|^2 prod_j (dx_j/2) for nu a multi-index in dimension nu

double* squared_L_sample_multi_d (size_t dim,mpa_multi_index* nu)
{
    double* x = uniform_sample_multi_d(dim);
    for (int j=0; j<nu->size; j++) {
        size_t i    = nu->data[j].index;
        size_t nu_i = nu->data[j].coordinate;
        x[i] = squared_L_sample_1d(nu_i);
    }
    return x;
}


#pragma mark -
#pragma mark sampling from 1d and multi_d squared Chebechev polynomial based measure

double squared_T_pdf_1d (size_t k, double x)
{
    double Tk_x = univariate_Chebychev_polynomial(k, x, NULL);
    if (k==0)
        return arcsine_pdf_1d(x);
    else
        return 2*Tk_x*Tk_x*arcsine_pdf_1d(x);
}


double squared_T_pdf_multi_d (size_t dim,mpa_multi_index* nu, double* x)
{
    double Tnu_x = multivariate_Chebychev_polynomial(nu, x, NULL);
    return pow(2.0, nu->size)*Tnu_x*Tnu_x*arcsine_pdf_multi_d(dim, x);
}

// sample from density |sqrt{2} T_k(x)|^2 dx/2 for k>=1
// using rejection method based on Berstein type inequality
// for chebyshev polynomials |T_k(x)| (1-x^2)^{1/4} <=1
double squared_T_sample_1d (size_t k)
{
    double c = 2;
    double x, u, fx, gx;
    int i=0;
    //printf("---------------------\n");
    while (1) {
        u  = 1.0*rand()/RAND_MAX;
        x =  arcsine_sample_1d();
        fx = squared_T_pdf_1d(k, x);
        gx = arcsine_pdf_1d(x);
        if (fx > c *gx*u)
            break;
        else{
            //     printf("------%d\n",i);
            i+=1;
        }
    }
    return x;
}

// sample from density 2^{|nu|_0} |T_nu(x)|^2 prod_j (dx_j/2) for nu a multi-index in dimension nu

double* squared_T_sample_multi_d (size_t dim,mpa_multi_index* nu)
{
    double* x = arcsine_sample_multi_d(dim);
    for (int j=0; j<nu->size; j++) {
        size_t i    = nu->data[j].index;
        size_t nu_i = nu->data[j].coordinate;
        x[i] = squared_T_sample_1d(nu_i);
    }
    return x;
}


#pragma mark -
#pragma mark sampling from 1d and multi_d christoffel Legendre realted measure over [-1,1]


double christoffel_Legendre_kernel_1d(size_t n, double x)
{
    gsl_vector* vec_Lx = gsl_vector_calloc(n);
    gsl_vector_set(vec_Lx, 0, 1);
    gsl_vector_set(vec_Lx, 1, x);
    for (int k=1; k<n-1; k++) {
        double y = (2.*k+1.)*x* gsl_vector_get(vec_Lx, k) - k* gsl_vector_get(vec_Lx, k-1);
        y /= k+1.;
        gsl_vector_set(vec_Lx, k+1, y);
    }
    
    for (int k=0; k<n; k++) {
        double Pk_x = gsl_vector_get(vec_Lx, k);
        gsl_vector_set(vec_Lx, k, sqrt(2*k+1)*Pk_x);
    }
    double ans=0.0;
    gsl_blas_ddot(vec_Lx, vec_Lx, &ans);
    return ans;
}

double christoffel_Legendre_pdf_1d (size_t n, double x)
{
    double ans = christoffel_Legendre_kernel_1d(n, x)/n;
    return ans/2;
    
}

/*
 with acceptance-rejection method using arcsine distribution as the proposal distribution
 We have verified that f(x) < c g(x) where c=(e/2)~1.35 for any n>=3 where
 f(x) = christoffel_Legendre_pdf_1d(n,x) and g(x) = arcsine_pdf_1d(x)
 For n=1 and n=2, we simply use the constant c = pi/2.
*/
double* iid_christoffel_Legendre_sample_1d (size_t n, size_t N)
{
    double* v = (double*) calloc(N, sizeof(double));
    srand( (unsigned int)time( NULL ) );
    for (int i=0; i<N; i++){
        double x , u, fx, gx;
        double c = (n>=3)? 0.5*M_E:0.5*M_PI;
        while (1) {
            u  = 1.0*rand()/RAND_MAX;
            x  = cos(3.14*rand()/RAND_MAX);
            fx = christoffel_Legendre_pdf_1d(n, x);
            gx = arcsine_pdf_1d(x);
            if (fx > c *gx*u){
                v[i] = x;
                break;
            }
        }
    }
    return v;
}

void set_Legendre_vector (mpa_monotone_graph* Lambda, size_t activeSize, double *y, gsl_vector* L)
{
    for (int i=0; i<activeSize; i++) {
        mpa_multi_index* nu = Lambda->nodes[i]->miPtr;
        double Ly = multivariate_Legendre_polynomial_L(nu, y, NULL);
        gsl_vector_set(L, i, Ly);
    }
}

double christoffel_Legendre_kernel_multi_d (size_t dim, mpa_monotone_graph* Lambda, size_t activeSize, double* y)
{
    gsl_vector* vec_Ly = gsl_vector_alloc(activeSize);
    set_Legendre_vector (Lambda, activeSize, y, vec_Ly);
    double ans=0.0;
    gsl_blas_ddot(vec_Ly, vec_Ly, &ans);
    return ans;
}

double christoffel_Legendre_pdf_multi_d (size_t dim, mpa_monotone_graph* Lambda, size_t activeSize, double* y)
{
    double ker_y = christoffel_Legendre_kernel_multi_d (dim, Lambda, activeSize, y);
    return ker_y/(activeSize*pow(2., dim));
}




#pragma mark -
#pragma mark sampling from skew christoffel Legendre realted measure over [-1,1]


double skew_christoffel_Legendre_kernel_multi_d (size_t dim, mpa_monotone_graph* Lambda, size_t activeSize, gsl_matrix* M, double* y)
{
    gsl_vector* vec_Ly = gsl_vector_alloc(activeSize);
    set_Legendre_vector (Lambda, activeSize, y, vec_Ly);
    return Quadratic_form(M, vec_Ly, vec_Ly);
}

double skew_christoffel_Legendre_pdf_multi_d (size_t dim, mpa_monotone_graph* Lambda, size_t activeSize, gsl_matrix* M, double* y)
{
    double LZL  = skew_christoffel_Legendre_kernel_multi_d (dim, Lambda, activeSize, M, y);
    double trace= matrix_trace(M);
    return LZL/(pow(2.,dim)*trace);
}


double* christoffel_Legendre_sample_multi_d(size_t dim, mpa_monotone_graph* Lambda, size_t activeSize)
{
    int idx = rand()%activeSize;
    return squared_L_sample_multi_d(dim, Lambda->nodes[idx]->miPtr);
}


double** iid_christoffel_Legendre_sample_multi_d (size_t dim, mpa_monotone_graph* Lambda, size_t activeSize, size_t N)
{
    double** v = (double**) calloc(N, sizeof(double*));
    for (int i=0; i<N; i++)
        v[i] = christoffel_Legendre_sample_multi_d( dim, Lambda, activeSize);
    return v;
}

/*
 sample from the density L(x) M L(x)/ tr(M)) where M is SDP. The matrix M is given
*/

double* skew_christoffel_Legendre_sample_multi_d(size_t dim, mpa_monotone_graph* Lambda, size_t activeSize, gsl_matrix* M)
{
    double* x=NULL;
    doublePair minmax = sym_matrix_minmax_eigenvalue(M);
    double lam_max = minmax.second;
    double trace=matrix_trace (M);
    double c = activeSize * lam_max / trace;
    
    double u, fx, gx;
    while (1) {
        u  = 1.0*rand()/RAND_MAX;
        x  = christoffel_Legendre_sample_multi_d  (dim, Lambda, activeSize);
        fx = skew_christoffel_Legendre_pdf_multi_d(dim, Lambda, activeSize, M, x);
        gx = christoffel_Legendre_pdf_multi_d     (dim, Lambda, activeSize, x);
        if (fx > c *gx*u){
            break;
        }
        else{
            num_rejection+=1;
            free(x);
        }
    }
    return x;
}




#pragma mark -
#pragma mark Algorithms for geneating samples/weights/gram_matrices to be used in weighed least squares


ls_datas* w_uniform_sample_multi_d (size_t dim, mpa_monotone_graph* Lambda, size_t activeSize, size_t N)
{
    gsl_matrix* G = gsl_matrix_calloc(activeSize, activeSize);
    gsl_vector* L = gsl_vector_calloc(activeSize);
    
    double** x = (double**) calloc(N, sizeof(double*));
    double* s  = (double*)  calloc(N, sizeof(double));
    gsl_matrix* D = gsl_matrix_calloc(N, activeSize);
    
    //srand( (unsigned int)time( NULL ) );
    for (int i=0; i<N; i++){
        x[i] = uniform_sample_multi_d(dim);
        set_Legendre_vector(Lambda, activeSize, x[i], L);
        gsl_matrix_set_row(D, i, L);
        s[i] = 1.0;
        matrix_add_outer_product (G, L, s[i]);
    }
    
    ls_datas* ans = ls_datas_alloc(dim, N);
    ls_datas_set_points (ans, x);
    ls_datas_set_weights(ans, s);
    ls_datas_set_designMatrix(ans, D);
    ls_datas_set_gramMatrix  (ans, G);
    return ans;
}


ls_datas* w_arcsine_sample_multi_d (size_t dim, mpa_monotone_graph* Lambda, size_t activeSize, size_t N)
{
    gsl_matrix* G = gsl_matrix_calloc(activeSize, activeSize);
    gsl_vector* L = gsl_vector_calloc(activeSize);
    
    double** x = (double**) calloc(N, sizeof(double*));
    double* s  = (double*)  calloc(N, sizeof(double));
    gsl_matrix* D = gsl_matrix_calloc(N, activeSize);
    
    //srand( (unsigned int)time( NULL ) );
    for (int i=0; i<N; i++){
        x[i] = arcsine_sample_multi_d(dim);
        set_Legendre_vector(Lambda, activeSize, x[i], L);
        gsl_matrix_set_row(D, i, L);
        double inv_s = pow(2.0, dim) * arcsine_pdf_multi_d(dim, x[i]);
        s[i] = 1.0/inv_s;
        matrix_add_outer_product (G, L, s[i]);
    }
        
    ls_datas* ans = ls_datas_alloc(dim, N);
    ls_datas_set_points (ans, x);
    ls_datas_set_weights(ans, s);
    ls_datas_set_designMatrix(ans, D);
    ls_datas_set_gramMatrix  (ans, G);
    return ans;
}

ls_datas* vanilla_CL_sample_multi_d (size_t dim, mpa_monotone_graph* Lambda, size_t activeSize, size_t N)
{
    //double  reject_cst = christoffel_Legendre_rejection_cst (dim,  Lambda,  activeSize);
    gsl_matrix* G = gsl_matrix_calloc(activeSize, activeSize);
    gsl_vector* L = gsl_vector_calloc(activeSize);
    
    double** x = (double**) calloc(N, sizeof(double*));
    double* s  = (double*)  calloc(N, sizeof(double));
    gsl_matrix* D = gsl_matrix_calloc(N, activeSize);
    
    //srand( (unsigned int)time( NULL ) );
    for (int i=0; i<N; i++){
        x[i] = christoffel_Legendre_sample_multi_d(dim, Lambda, activeSize);
        set_Legendre_vector(Lambda, activeSize, x[i], L);
        gsl_matrix_set_row(D, i, L);
        double norm2_L = gsl_blas_dnrm2(L);
        s[i] = activeSize/pow(norm2_L,2);
        matrix_add_outer_product (G, L, s[i]);
    }
    
    ls_datas* ans = ls_datas_alloc(dim, N);
    ls_datas_set_points (ans, x);
    ls_datas_set_weights(ans, s);
    ls_datas_set_designMatrix(ans, D);
    ls_datas_set_gramMatrix  (ans, G);
    return ans;
}

#pragma mark -
#pragma mark Implementation of the effective resistance algorithms of the paper


typedef struct{
    double p0;
    double p1;
    double p2;
}doubleTriple;



/*
 Compute the quantities |L(x)|^2, L(x) invA L(x) and |invA L(x)|^2 where
 invA is the inverse of matrix A. The matrix A is given, not its inverse, and
 also A is not in plain format, but its cholesky decomposition provided by
 gsl_linalg_cholesky_decomp(M).
*/

doubleTriple three_special_squared_norms(const gsl_matrix* A, const gsl_vector* x)
{
    gsl_vector *y = gsl_vector_alloc(x->size);
    gsl_vector_memcpy(y, x);
    
    doubleTriple ans = {0.0,0.0,0.0};
    ans.p0  = dot_product(y, y);
    gsl_blas_dtrsv (CblasLower, CblasNoTrans, CblasNonUnit, A, y);
    ans.p1 = dot_product(y, y);
    gsl_blas_dtrsv (CblasLower,   CblasTrans, CblasNonUnit, A, y);
    ans.p2 = dot_product(y, y);
    gsl_vector_free(y);
    return ans;
}


/*
************************************************************************************
 Algo 3.1 variable barrier push (vbp) sampling
************************************************************************************
*/


/*
 sample from the density L (M^{-1} + gamma Id/n ) L/ (tr(M^{-1}) + gamma)
 where M is SDP. The matrix M is given, not its inverse, and also M is not in plain
 format, but its cholesky decomposition provided by gsl_linalg_cholesky_decomp(M).
 Also are inputs the largest eigenvalue of M^{-1}
 */

typedef struct{
    double *x;
    gsl_vector* Lx;
    doubleTriple norms;
} point_and_norms;


point_and_norms skew_CL_sample_multi_d(size_t dim, mpa_monotone_graph* Lambda, size_t activeSize,
                                       gsl_matrix* M, double lam_max_inv_M, double gamma)
{
    double* x=NULL;
    double lam_max = lam_max_inv_M + gamma/activeSize;
    double c = lam_max ;
    gsl_vector* Lx = gsl_vector_alloc(activeSize);
    doubleTriple dt;
    
    
    double u, fx, gx;
    while (1) {
        u  = 1.0*rand()/RAND_MAX;
        x  = christoffel_Legendre_sample_multi_d(dim, Lambda, activeSize);
        set_Legendre_vector(Lambda, activeSize, x, Lx);
        dt =  three_special_squared_norms(M, Lx);
        
        fx = dt.p1 + (gamma/activeSize)*dt.p0; // not pdf
        gx = dt.p0; // not pdf
        if (fx > c *gx*u){
            break;
        }
        else{
            num_rejection+=1;
            free(x);
        }
    }
    //printf("%lu, ",num_rejection);
    //num_rejection=0;

    point_and_norms ans = {x,Lx,dt};
    return ans;
}


ls_datas* vbp_CL_sample_multi_d(size_t dim, mpa_monotone_graph* Lambda, size_t activeSize,
                                double eps, double gamma, size_t N)
{
        
    double** x = (double**) calloc(N, sizeof(double*));
    double* s  = (double*)  calloc(N, sizeof(double));
    gsl_matrix* A  = gsl_matrix_calloc(activeSize, activeSize);
    gsl_matrix* D  = gsl_matrix_calloc(N, activeSize);
    gsl_matrix* AA = gsl_matrix_calloc(activeSize, activeSize);
        
    double l, trY,trZ, delta;
    double eta = eps/(1.-eps);
    
    l = - 1.0*activeSize;
    trY = 1.0;
        
    for (int i=0; i<N; i++){
        delta = eps/(trY+gamma);
        l+=delta;
        
        gsl_matrix_memcpy(AA, A); // we copy A, and work on it
        matrix_add_multiple_eye(AA, -l);
        gsl_vector* eigenAA = sym_matrix_eigenvalues (AA);
        trZ=0.0;
        for (int i=0; i<eigenAA->size; i++)
            trZ += 1./gsl_vector_get(eigenAA, i);

        double lam_min = gsl_vector_min(eigenAA);
        double lam_max_Z = 1.0/lam_min;
        gsl_vector_free(eigenAA);
        
        gsl_linalg_cholesky_decomp(AA);
        point_and_norms pn = skew_CL_sample_multi_d(dim, Lambda, activeSize,
                                                    AA, lam_max_Z, gamma);
            
        x[i]=pn.x;
        gsl_matrix_set_row(D, i, pn.Lx);
        
        double LL   = pn.norms.p0;
        double LZL  = pn.norms.p1;
        double LZZL = pn.norms.p2;
        
        double rho_x  = LZL + (gamma/activeSize) * LL;
        s[i] = eta/rho_x;
        matrix_add_outer_product (A, pn.Lx, s[i]);
            
        trY = trZ - s[i]*LZZL/ (1+s[i]*LZL);
    }
        
    gsl_matrix_free(AA);
    
    ls_datas* ans = ls_datas_alloc(dim, N);
    ls_datas_set_points (ans, x);
    ls_datas_set_weights(ans, s);
    ls_datas_set_designMatrix(ans, D);
    ls_datas_set_gramMatrix  (ans, A);
    return ans;
}

/*
************************************************************************************
 Algo 4.1
************************************************************************************
*/


/*
sample from the probability density associated with non-normalized measure w(x) 1_{w(x)>=treshold}
where w(x) = L(x) (Z^2/diff_trace -Z) L(x) and Z = M^{-1}, M is SDP.
The matrix M is given, not its inverse, and also M is not in plain format, but its
cholesky decomposition provided by gsl_linalg_cholesky_decomp(M).
Also are inputs the largest eigenvalue of M^{-1}
 */


point_and_norms skew_CL_sample_multi_d_conditionned(size_t dim, mpa_monotone_graph* Lambda, size_t activeSize,
                                                    gsl_matrix* M, double diff_trace, double lam_max_Z, double treshold)
{
    double* x=NULL;
    double lam_max_W = lam_max_Z*lam_max_Z/diff_trace - lam_max_Z;// we know it is >0
    double c = lam_max_W ; // lambda max of matrix  W = Z*Z/diff_trace - Z
    gsl_vector* Lx = gsl_vector_alloc(activeSize);
    doubleTriple dt;
    
    
    double u, fx, gx;
    while (1) {
        u  = 1.0*rand()/RAND_MAX;
        x  = christoffel_Legendre_sample_multi_d(dim, Lambda, activeSize);
        set_Legendre_vector(Lambda, activeSize, x, Lx);
        dt =  three_special_squared_norms(M, Lx);
        
        fx = dt.p2/diff_trace - dt.p1; // not a value of pdf
        gx = dt.p0; // not a value of pdf
        if (fx > treshold && fx > c *gx*u){
            break;
        }
        else{
            num_rejection+=1;
            free(x);
        }
    }
    //printf("%lu, ",num_rejection);
    //num_rejection=0;

    point_and_norms ans = {x,Lx,dt};
    return ans;
}


/*
************************************************************************************
 Algo 4.1 fixed barrier push (fbp) sampling
************************************************************************************
*/


ls_datas* fbp_CL_sample_multi_d(size_t dim, mpa_monotone_graph* Lambda, size_t activeSize,
                                double delta, double kappa, size_t N)
{
    
    
    double** x = (double**) calloc(N, sizeof(double*));
    double* s  = (double*)  calloc(N, sizeof(double));
    gsl_matrix* A = gsl_matrix_calloc(activeSize, activeSize);
    gsl_matrix* D = gsl_matrix_calloc(N, activeSize);
    gsl_matrix* AA = gsl_matrix_calloc(activeSize, activeSize);
    
    double threshold = kappa*(1.0-delta)/delta;

    double l,trY,trZ ;
    gsl_matrix_set_all(A, 0.0);
    l = - 1.0*activeSize;
    trY = 1; // it will stays equal to 1
    for (int i=0; i<N; i++){
        l+=delta;
        gsl_matrix_memcpy(AA, A);
        matrix_add_multiple_eye(AA, -l);
        gsl_vector* eigenAA = sym_matrix_eigenvalues (AA);
        double lam_min = gsl_vector_min(eigenAA);

        trZ=0.0;
        for (int i=0; i<eigenAA->size; i++)
            trZ += 1./gsl_vector_get(eigenAA, i);
        
        gsl_vector_free(eigenAA);
        
        double diff_trace = trZ-trY;
        double lam_max_Z = 1.0/lam_min;
        gsl_linalg_cholesky_decomp(AA);
        
        point_and_norms pn = skew_CL_sample_multi_d_conditionned(dim, Lambda, activeSize,
                                                                 AA, diff_trace, lam_max_Z, threshold);


        x[i]=pn.x;
        gsl_matrix_set_row(D, i, pn.Lx);
        double w_x  = pn.norms.p2/diff_trace - pn.norms.p1;
        s[i] = 1.0/w_x;
        matrix_add_outer_product (A, pn.Lx, s[i]);
            
        double LZL  = pn.norms.p1;
        double LZZL = pn.norms.p2;
        trY = trZ - s[i]*LZZL/ (1+s[i]*LZL); // not needed, just to make sure it stays equal to 1
    }
                
    gsl_matrix_free(AA);
    ls_datas* ans = ls_datas_alloc(dim, N);
    ls_datas_set_points (ans, x);
    ls_datas_set_weights(ans, s);
    ls_datas_set_designMatrix(ans, D);
    ls_datas_set_gramMatrix  (ans, A);
    return ans;
}
