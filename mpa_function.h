/*
 *  mpa_function.h
 *  mpa
 *
 *  Created by Abdellah CHKIFA on 07/03/11.
 *  Copyright 2011 Oliva Electronics. All rights reserved.
 *
 */


#ifndef mpa_function_solver_h
#define mpa_function_solver_h

typedef struct
{
    //size_t dim;
	void* data;
	double (*evaluate)(void* arg,void* data);
}mpa_function;

typedef struct
{
    size_t dim;
    void* datas;
    int (*evaluate)(const double* input, size_t dim, void* datas, void* output);
} query_solver;



#endif /* mpa_function_solver_h */
