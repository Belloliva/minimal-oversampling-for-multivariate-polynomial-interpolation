/*
 *  mpa_combinatorics.h
 *  mpa
 *
 *  Created by Abdellah CHKIFA on 02/03/11.
 *  Copyright 2011 oliva Electronics. All rights reserved.
 *
 */

#include <stdlib.h>

#ifndef __mpa_combinatorics_h__
#define __mpa_combinatorics_h__

/*!
 * \brief an integer factorial.
 
 * Return the factorial of n, i.e. n!=n*(n-1)*...*2
 * @param n the integer
 * @return the factorial value
 */
size_t mpa_factorial(size_t n);

/*!
 * \brief binomial coefficient.
 
 * Return the binomial coefficient of n and k. \f$ {n \choose k} =\frac{n!}{k!(n-k)!}\f$
 * @param n the first integer.
 * @param k the second integer.
 * @return the binomial coefficient value.
 */
size_t mpa_binomial_coefficient(size_t n, size_t k);

#endif /*__mpa_combinatorics_h__*/
