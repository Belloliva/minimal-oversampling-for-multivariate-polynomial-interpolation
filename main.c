//
//  main_scalar_functions.c
//  ff_tests
//
//  Created by abdellah chkifa on 11/27/23.
//


#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#include "mpa_multi_index.h"
#include "mpa_function.h"
#include "mpa_linalg.h"
#include "mpa_sampling.h"
#include "Least_squares_solver.h"

size_t num_rejection=0;


void print_vector_terminal_d (double* vec, size_t num_elements)
{
    printf("[");
    for(int j=0; j<num_elements-1; j++)
        printf("%.6f,",vec[j]);
    printf("%.6f]; \n\n",vec[num_elements-1]);
}

void print_vector_terminal_lu (size_t* vec, size_t num_elements)
{
    printf("[");
    for(int j=0; j<num_elements-1; j++)
        printf("%lu,",vec[j]);
    printf("%lu]; \n\n",vec[num_elements-1]);
}


int multivariate_generating_function_Legendre(const double* input, size_t dim, void* data, void* output){
    const double* x = (const double* )input;
    double* fx = (double*)output;
    double* rho = (double*)data;
    
    double prod =1.;
    for(int j=0;j<dim;j++){
        double xj = x[j];
        double t  = rho[j];
        prod *= 1./sqrt(1- 2*t*xj + t*t) ;
    }
    *fx =  prod ;//+0.01*(2.*rand()/RAND_MAX-1.);
    return EXIT_SUCCESS;
}



double multivariate_generating_function_coefficient_P (void* arg,void* data)
{
    mpa_multi_index* nu = (mpa_multi_index*)arg;
    double* rho = (double*)data;
    
    double prod = 1.;
    for (int j=0; j<nu->size; j++){
        size_t ix_j = nu->data[j].index;
        size_t nu_j = nu->data[j].coordinate;
        prod *= pow(rho[ix_j],nu_j);
    }
    return prod;
}


double multivariate_generating_function_coefficient_L (void* arg,void* data)
{
    mpa_multi_index* nu = (mpa_multi_index*)arg;
    double* rho = (double*)data;
    
    double prod = 1;
    for (int j=0; j<nu->size; j++){
        size_t ix_j = nu->data[j].index;
        size_t nu_j = nu->data[j].coordinate;
        prod *= pow(rho[ix_j],nu_j);
        prod /= sqrt(2.0*nu_j+1.0);
    }
    return prod;
}


int experiment1 (void){
    time_t start = time(NULL);
    
    size_t dim = 4;
    double* rho = (double*)calloc(dim, sizeof(double));
    for(int j=0;j<dim;j++)
        rho[j] = 0.9 - 0.1*j;
    
    // First, we define the solver used for querying the target function
    query_solver qsolver = {dim, rho, multivariate_generating_function_Legendre};
    
    // then create the monotone graph (i.e. lower set) associated
    // with the n largest values of its Legendre coefficients
    
    size_t N_basis = 128;
    mpa_function func = {rho, multivariate_generating_function_coefficient_L};
    mpa_monotone_graph* Lambda = mpa_monotone_graph_envellope(dim, N_basis, func);
    mpa_monotone_graph_resize(Lambda, N_basis);
    
    for (int i=0; i<Lambda->size; i++) {
        double* d = (double*)Lambda->nodes[i]->info;
        free(d);
        Lambda->nodes[i]->info=NULL;
    }
    
    // Then we compute the best L2 norm associated with approximating the function
    // in the polynomial space associated with the lower set Lambda.
    
    double value_at_1 =1;
    double L2_norm =1;
    
    for (int j=0; j<dim; j++){
        double x = rho[j];
        L2_norm   *= sqrt((log(1+x)-log(1-x))/(2*x));
        value_at_1 *= 1.0/(1-x);
    }
    
    double* exact_coef_P = (double*) calloc(N_basis, sizeof(double));
    double* exact_coef_L = (double*) calloc(N_basis, sizeof(double));
    
    double value_at_1_Lambda = 0.0;
    double L2_norm_squared_Lambda = 0.0;
    
    for (int i=0; i<N_basis; i++) {
        mg_node* node = Lambda->nodes[i];
        mpa_multi_index* nu = node ->miPtr;
        
        exact_coef_P[i] = multivariate_generating_function_coefficient_P(nu, rho);
        exact_coef_L[i] = multivariate_generating_function_coefficient_L(nu, rho);
        
        value_at_1_Lambda       += exact_coef_P[i];
        L2_norm_squared_Lambda  += pow(exact_coef_L[i], 2.0);
    }
    
    double best_L2_error_squared = pow(L2_norm,2.0) - L2_norm_squared_Lambda;
    double value_at_1_error = value_at_1 - value_at_1_Lambda; // for approximating at
    
    
    printf("L2 error of approxiation  =%f  \n",sqrt(best_L2_error_squared) );
    printf("error of approximation at 1 =%f \n",value_at_1_error );
    
    /*
     ----------------------------------------------------------------------------
     Below is the experiment, we simply uncomment a line ls_datas* ....=
     and rerun the code to get the outputs for the specified sampling strategy
     ----------------------------------------------------------------------------
     */
    
    size_t numPoints = floor(2.0*N_basis);
    size_t num_trials = 400;
    
    srand( (unsigned int)time( NULL ) );
    
    double* L2_error        = (double*) calloc(num_trials, sizeof(double));
    double* cond_number     = (double*) calloc(num_trials, sizeof(double));
    size_t* num_rejections  = (size_t*) calloc(num_trials, sizeof(size_t));
    
    // the value of r for algorithm 3.1, and associated epsilon and gamma
    double r1 = (numPoints+1.0)/N_basis;
    double eps = pow(r1, -0.25);
    double gamma = pow(r1, 0.5) - pow(r1, 0.25);
    
    // the value of r for algorithm 4.1, and associated delta and kappa
    double r2 = numPoints/(N_basis-1.0);
    double delta = 1/sqrt(r2);
    double kappa = 0.5;

    printf("\n\n Simulation started\n\n");

    for (int k=0; k<num_trials; k++) {
        if (k%10==0)
            printf("trials 1 to %d done \n",k);
        
        /* iid sampling */
        //ls_datas* ls_datas0 = w_uniform_sample_multi_d(dim, Lambda, N_basis, numPoints);
        //ls_datas* ls_datas0 = w_arcsine_sample_multi_d(dim, Lambda, N_basis, numPoints);
        ls_datas* ls_datas0 = vanilla_CL_sample_multi_d (dim, Lambda, N_basis, numPoints);
        
        /*
         history dependant sampling proposed by the paper
         vbp_CL_sample_multi_d is algorithm 3.1 of the paper
         fbp_CL_sample_multi_d is algorithm 4.1 of the paper
         */
        
        //ls_datas* ls_datas0 = vbp_CL_sample_multi_d(dim, Lambda, N_basis, eps, gamma, numPoints);
        //ls_datas* ls_datas0 = fbp_CL_sample_multi_d(dim, Lambda, N_basis, delta, kappa, numPoints);
        
        num_rejections[k] = num_rejection;
        num_rejection = 0;
        cond_number[k] = sym_matrix_condition_number(ls_datas0->gramMatrix);
        
        set_ls_solution(ls_datas0, Lambda, N_basis, qsolver);
        ls_datas_free(ls_datas0);
        
        double L2_error_squared = 0;
        for (int i=0; i<N_basis; i++) {
            mg_node* node = Lambda->nodes[i];
            double* info   = (double*)node->info;
            double ls_coef = *info;
            L2_error_squared += pow(ls_coef - exact_coef_L[i], 2.0);
        }
        L2_error[k] = sqrt(L2_error_squared + best_L2_error_squared)/sqrt(best_L2_error_squared);
    }
    /*
    ------------------------------------------------------------------------------------------------
    ------------------------------------------------------------------------------------------------
    */
    
    printf("condition_number=");
    print_vector_terminal_d(cond_number,num_trials);
    
    printf("rejection_number=");
    print_vector_terminal_lu(num_rejections,num_trials);
    
    printf("L2_relative_error=");
    print_vector_terminal_d(L2_error,num_trials);
    
    //printf("number of rejection is =%lu \n ", num_rejection );
    printf("Time taken by the runs is %f secondes\n",difftime(time(NULL),start));
    return 0;
}

int experiment2 (void){
    time_t start = time(NULL);
    
    size_t dim = 4;
    double* rho = (double*)calloc(dim, sizeof(double));
    for(int j=0;j<dim;j++)
        rho[j] = 0.9 - 0.1*j;
    
    // First, we define the solver used for querying the target function
    query_solver qsolver = {dim, rho, multivariate_generating_function_Legendre};
    
    // then create the monotone graph (i.e. lower set) associated
    // with the n largest values of its Legendre coefficients
    
    size_t N_basis = 128;
    mpa_function func = {rho, multivariate_generating_function_coefficient_L};
    mpa_monotone_graph* Lambda = mpa_monotone_graph_envellope(dim, N_basis, func);
    mpa_monotone_graph_resize(Lambda, N_basis);
    
    for (int i=0; i<Lambda->size; i++) {
        double* d = (double*)Lambda->nodes[i]->info;
        free(d);
        Lambda->nodes[i]->info=NULL;
    }
    
    // Then we compute the best L2 norm associated with approximating the function
    // in the polynomial space associated with the lower set Lambda.
    
    double value_at_1 =1;
    double L2_norm =1;
    
    for (int j=0; j<dim; j++){
        double x = rho[j];
        L2_norm   *= sqrt((log(1+x)-log(1-x))/(2*x));
        value_at_1 *= 1.0/(1-x);
    }
    
    double* exact_coef_P = (double*) calloc(N_basis, sizeof(double));
    double* exact_coef_L = (double*) calloc(N_basis, sizeof(double));
    
    double value_at_1_Lambda = 0.0;
    double L2_norm_squared_Lambda = 0.0;
    
    for (int i=0; i<N_basis; i++) {
        mg_node* node = Lambda->nodes[i];
        mpa_multi_index* nu = node ->miPtr;
        
        exact_coef_P[i] = multivariate_generating_function_coefficient_P(nu, rho);
        exact_coef_L[i] = multivariate_generating_function_coefficient_L(nu, rho);
        
        value_at_1_Lambda       += exact_coef_P[i];
        L2_norm_squared_Lambda  += pow(exact_coef_L[i], 2.0);
    }
    
    double best_L2_error_squared = pow(L2_norm,2.0) - L2_norm_squared_Lambda;
    double value_at_1_error = value_at_1 - value_at_1_Lambda; // for approximating at
    
    
    printf("L2 error of approxiation  =%f  \n",sqrt(best_L2_error_squared) );
    printf("error of approximation at 1 =%f \n",value_at_1_error );
    
    /*
     ------------------------------------------------------------------------------------------------
     Below is the experiment, we simply uncomment a line ls_datas* ....=
     and run the code
     ------------------------------------------------------------------------------------------------
     */
    
    size_t numPoints = N_basis;
    
    size_t N_p = 40;
    size_t num_trials = 100;
    printf("\n\n L2_relative_error = np.zeros((%lu,%lu))\n\n",N_p,num_trials);


    srand( (unsigned int)time( NULL ) );
    for (int p=0; p<N_p; p++) {
        numPoints = N_basis + 2*p;
        
        double* L2_error        = (double*) calloc(num_trials, sizeof(double));
        //double* cond_number     = (double*) calloc(num_trials, sizeof(double));
        //size_t* num_rejections  = (size_t*) calloc(num_trials, sizeof(size_t));
        
        // the value of r for algorithm 3.1, and epsilon and gamma
        double r1 = (numPoints+1.0)/N_basis;
        double eps = pow(r1, -0.25);
        double gamma = pow(r1, 0.5) - pow(r1, 0.25);
        
        // the value of r for algorithm 4.1, and delta and kappa
        double r2 = numPoints/(N_basis-1.0);
        double delta = 1/sqrt(r2);
        double kappa = 0.5;
        
        for (int k=0; k<num_trials; k++) {
            
            /* iid sampling */
            //ls_datas* ls_datas0 = w_uniform_sample_multi_d(dim, Lambda, N_basis, numPoints);
            //ls_datas* ls_datas0 = w_arcsine_sample_multi_d(dim, Lambda, N_basis, numPoints);
            ls_datas* ls_datas0 = vanilla_CL_sample_multi_d (dim, Lambda, N_basis, numPoints);
            
            /*
             history dependant sampling proposed by the paper
             vbp_CL_sample_multi_d is algorithm 3.1 of the paper
             fbp_CL_sample_multi_d is algorithm 4.1 of the paper
             */
            
            //ls_datas* ls_datas0 = vbp_CL_sample_multi_d(dim, Lambda, N_basis, eps, gamma, numPoints);
            //ls_datas* ls_datas0 = fbp_CL_sample_multi_d(dim, Lambda, N_basis, delta, kappa, numPoints);
            
            /*
            num_rejections[k] = num_rejection;
            num_rejection = 0;
            cond_number[k] = sym_matrix_condition_number(ls_datas0->gramMatrix);
            */
            
            set_ls_solution(ls_datas0, Lambda, N_basis, qsolver);
            ls_datas_free(ls_datas0);
            
            double L2_error_squared = 0;
            for (int i=0; i<N_basis; i++) {
                mg_node* node = Lambda->nodes[i];
                double* info   = (double*)node->info;
                double ls_coef = *info;
                L2_error_squared += pow(ls_coef - exact_coef_L[i], 2.0);
            }
            L2_error[k] = sqrt(L2_error_squared + best_L2_error_squared)/sqrt(best_L2_error_squared);
            
        }
        /*
         ------------------------------------------------------------------------------------------------
         We then print for every p the error values, the print is to be used by python to produce figures
         ------------------------------------------------------------------------------------------------
         */
        
        /*
        printf("condition_number[%d,:]=",p);
        print_vector_terminal_d(cond_number,num_trials);
        
        printf("rejection_number[%d,:]=",p);
        print_vector_terminal_lu(num_rejections,num_trials);
        */
        
        printf("L2_relative_error[%d,:]=",p);
        print_vector_terminal_d(L2_error,num_trials);
    }
    printf("Time taken by the runs is %f secondes\n",difftime(time(NULL),start));
    return 0;
}



int main (int argc, const char * argv[]){
    experiment1 ();
    //experiment2 ();
    return 0;
}
