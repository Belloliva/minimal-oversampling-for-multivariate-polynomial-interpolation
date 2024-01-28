/*
 *  mpa_combinatorics.c
 *  bulk_chase
 *
 *  Created by Abdellah CHKIFA on 02/03/11.
 *  Copyright 2011 oliva Electronics. All rights reserved.
 *
 */


/* multiple functions in combinatorics*/
#include "mpa_combinatorics.h"
#include <assert.h>


// factorial for small values of n
size_t mpa_factorial(size_t n)
{
    return (n==0)?1:n* mpa_factorial(n-1);
}

size_t mpa_binomial_coefficient(size_t n,size_t k)
{
    if (k==0)
        return 1;
    if (n==0)
        return 0;
    int res=1;
    int i=0;
    while(i++<k)
        res*=(n-k+i);
    
    return res/mpa_factorial(k);
    // We can also use Pascal triangle, but has to be well coded to prevent slow recursive call
    // return mpa_binomial_coefficient(n-1,k-1)+mpa_binomial_coefficient(n-1,k);
    // The best solution is to cache the triangle for n in {0,...,N} and k in {0,...,K}
    // for N and K big enough, then simply use the table.
}
