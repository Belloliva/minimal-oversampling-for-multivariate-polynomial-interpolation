/*
 *  mpa_multi_index.h
 *  mpa
 *
 *  Created by Abdellah CHKIFA on 17/02/11.
 *  Copyright 2011 oliva Electronics. All rights reserved.
 *
 */

#ifndef __mpa_multi_index_h__
#define __mpa_multi_index_h__

#include "mpa_stdlib.h"
#include <stdio.h>
#include <math.h>

/*! \struct mpa_pairIC
 \brief pair (index,coordinate).
 
 A pair of two unsigned int used to construct the multi-index struct.
 */
typedef struct
{
    size_t index;/*!< The coordinate index. */
	size_t coordinate;/*!< The coordinate value.*/
}
mpa_pairIC;//a pair (index,coordinate), to represent active coordinates

/*! \struct mpa_multi_index
 \brief multi-index notation.
 
 The mathematical notation of multi-indices simplifies formulae used in multivariable calculus, partial differential equations and the theory of distributions,
 by generalising the concept of an integer index to an ordered tuple of indices. see http://en.wikipedia.org/wiki/Multi-index_notation.
 for a multi-index \f$\nu=(1,0,3,5,0)\f$, the length is 5, the size is 3 and the data is {(0,1),(2,3),(3,5)}.
 */
typedef struct
{
    size_t size;/*!<The number of the multi-index active coordinates.*/
    size_t length;/*!<The length of the multi_index.*/
	mpa_pairIC* data;/*!<The array of the multi-index active coordinates, the index-coordinate pairs are sorted in an ascending order comparing the indices, and coordinates are never equal to zero.*/
}
mpa_multi_index;

#pragma mark -
#pragma mark memroy management

/*!
 * Creates a multi-index, returning a pointer to a newly initialized multi-index struct. A new data is allocated for the index-coordinate pairs of the multi-index, and
 * stored in the data component of the multi-index struct. The data is “owned” by the multi-index, and will be deallocated when the multi-index is deallocated.
 
 * @param n the multi-index size.
 * @param l the multi-index length.
 * @return pointer the new allocated multi-index
 */
mpa_multi_index *mpa_multi_index_alloc (const size_t n,const size_t l);
/*!
 * Free a previously allocated multi-index mi.
 * @param mi the multi-index to be freed.
 */
void mpa_multi_index_free(mpa_multi_index* mi);
/*!
 Copy the elements of the multi-index src into the multi-index dest. The two multi-indices must have the same length and the same size.
 * @param src the multi-index to copy from.
 * @param dest the multi-index to  copy to.
*/ 
void mpa_multi_index_memcpy (mpa_multi_index * dest,const mpa_multi_index * src);

#pragma mark -
#pragma mark getters
/*!
Return a pointer to the i-th index-coordinate pair of the multi-index mi. A null pointer is returned if i is not an active coordinate of mi
 *@param mi the multi-index.
 *@param i the index.
 *@return pointer to the i-th index-coordinate pair mi.
*/
mpa_pairIC* mpa_multi_index_get(const mpa_multi_index * mi,size_t i);

/*!
 *Return a pointer to a newly allocated array of length the multi-index mi length and holding mi coordinates in positions where mi is active and 0
 *elsewhere. For a multi-index of length 5, of size 3 and the data {(0,1),(2,3),(3,5)}. then the result is the array {1,0,3,5,0}.
 * @param mi the multi-index .
 * @return pointer the new allocated array
 */
size_t* mpa_multi_index_vector(const mpa_multi_index * mi);

#pragma mark -
#pragma mark simple functions
/*!
 *Return the absolute value of the multi-index, the sum of all it's coordinates. For a multi-index  \f$\nu=(\nu_0,\nu_1,...,\nu_{n-1})\f$, then the result is  \f$|\nu_0|+|\nu_1|+...+|\nu_{n-1}|\f$.
 * @param mi the multi-index.
 * @return the absolute value.
 */
size_t mpa_multi_index_abs(const mpa_multi_index * mi);

/*! 
 Return a weighted sum of the multi-index coordinates. The weight must be an array of the same length as the multi-index.
 For a multi-index  \f$\nu=(\nu_0,\nu_1,...,\nu_{n-1})\f$, and \f$weight=(w_0,...,w_{n-1})\f$ then the result is  \f$|\nu_0|w_0+|\nu_1|w_1+...+|\nu_{n-1}|w_{n-1}\f$.
 * @param mi the multi-index.
 * @param weight the weight to be used.
 * @return the weighted sum.
 */
double mpa_multi_index_dot_vector(const mpa_multi_index * mi,const double* weight);

/*!
 *Return the maximum of the multi-index coordinates. For a multi-index  \f$\nu=(\nu_0,\nu_1,...,\nu_{n-1})\f$, the result is \f$max_{0\leq i \leq n-1}(\nu_i)\f$.
 * @param mi the multi-index.
 * @return the maximum.
 */
size_t mpa_multi_index_max(const mpa_multi_index * mi);

/*!
 *Return the minimum of the multi-index coordinates. For a multi-index  \f$\nu=(\nu_0,\nu_1,...,\nu_{n-1})\f$, the result is \f$min_{0\leq i \leq n-1}(\nu_i)\f$.
 * @param mi the multi-index.
 * @return the minimum.
 */
size_t mpa_multi_index_min(const mpa_multi_index * mi);//min coordinate;

/*!
 Write the multi-index to the console. The multi-index coordinates are printed between two parentheses and separated by commas. The multi-index of length 5,
 size 3 and having data {(0,1),(2,3),(3,5)} is printed as \f$(1,0,3,5,0)\f$.
* @param mi the multi-index
*/ 
void mpa_multi_index_printf(const mpa_multi_index * mi);

/*!
 Write the multi-index to a stream. The multi-index coordinates are printed horizontally and separated by white spaces. The multi-index of length 5, size 3 and data={(0,1),(2,3),(3,5)} is printed as  \f$1~ 0~ 3~ 5~ 0~\f$.
 * @param mi the multi-index.
 * @param fp the stream.
 */ 
void mpa_multi_index_fprintf(FILE *fp,const mpa_multi_index* mi);

#pragma mark -
#pragma mark functions that give boolean values

/*!
Check if the multi-index is the zero multi-index. For \f$\nu=(\nu_0,\nu_1,...,\nu_{n-1})\f$, this mean that \f$\forall ~0\leq i \leq n-1~ \nu_i\= 0\f$. Return 1 if true and 0 if not.
 * @param mi the multi-index.
 * @return the boolean value
 */
int mpa_multi_index_isZero(const mpa_multi_index * mi);

/*!
Check if the multi-index is a unit vector, also called Kronecker vector. For \f$\nu=(\nu_0,\nu_1,...,\nu_{n-1})\f$, this mean that \f$\nu_j= 0\f$ for some \f$j\f$ and all the other coordinates are zeros. Return 1 if true and 0 if not.
 * @param mi the multi-index.
 * @return the boolean value
 */
int mpa_multi_index_isUnitVector(const mpa_multi_index * mi);

/*!
Check if the multi-index lies inside a hypercube. For \f$\nu=(\nu_0,\nu_1,...,\nu_{n-1})\f$, this mean that \f$\forall ~0\leq i \leq n-1~ \nu_i\leq side\f$. Return 1 if true and 0 if not.
 * @param mi the multi-index.
 * @param side the hypercude side.
 * @return the boolean value
 */ 
int mpa_multi_index_isInsideHypercube(const mpa_multi_index * mi,size_t side);

/*!
 Check if the multi-index lies inside a simplex of equation \f$ <w,\nu> \leq k\f$. Return 1 if true and 0 if not.
 * @param mi the multi-index.
 * @param w the weight used.
 * @param k the parameter k used.
 * @return the boolean value
 */ 
int mpa_multi_index_isInsideSimplex(const mpa_multi_index * mi,double* w, double k);

#pragma mark -
#pragma mark successor and predecessor
/*!
 Check if two multi-indices are equal. The two multi-indices must have the same length.  Return 1 if true and 0 if not.
 * @param mi1 the first multi-index.
 * @param mi2 the second multi-index.
 * @return the boolean value
 */
int mpa_multi_index_equal(const mpa_multi_index * mi1,const mpa_multi_index * mi2);


/*!
 Check if the multi-index is coordinates-wise less than a vector. The multi-index and vector must have the same length.
 If \f$\nu=(\nu_0,\nu_1,...,\nu_{n-1})\f$ and \f$ v=(v_0,v_1,...,v_{n-1})\f$, then
 \f$\nu<v\f$ if and only if \f$\forall 0\leq i \leq n-1 ~~ \nu_i\leq v_i\f$ and \f$\nu\neq v\f$.
 Return 1 if true and 0 if not.
 
 * @param mi the  multi-index.
 * @param vec the vector.
 * @return the boolean value
 */
int mpa_multi_index_less_vector (const mpa_multi_index * mi,const size_t* vec);

int mpa_multi_index_lessOrEqual (const mpa_multi_index * mi1,const mpa_multi_index * mi2);

int mpa_multi_index_less (const mpa_multi_index * mi1,const mpa_multi_index * mi2);

/*!
 Check if the multi-index is coordinates-wise less than the multi-index mi2. The two multi-indices must have the same length. If
 \f$\nu=(\nu_0,\nu_1,...,\nu_{n-1})\f$ and \f$\mu=(\mu_0,\mu_1,...,\mu_{n-1})\f$, then
 \f$\nu<\mu\f$ if and only if \f$\forall 0\leq i \leq n-1 ~~ \nu_i\leq \mu_i\f$ and
 \f$\nu\neq\mu\f$. Return 1 if true and 0 if not.
 
 * @param mi1 the first multi-index.
 * @param mi2 the second multi-index.
 * @return the boolean value
 */
int mpa_multi_index_less_helper (const mpa_multi_index * mi1,const mpa_multi_index * mi2,size_t* helper);

/*!
 Check if the multi-indicex is lexically less than the multi-index mi2. The two multi-indices must have the same length. If \f$\nu=(\nu_0,\nu_1,...,\nu_{n-1})\f$ and \f$\mu=(\mu_0,\mu_1,...,\mu_{n-1})\f$, then \f$\nu<_l\mu\f$ if and only if \f$\nu_0>\mu_0\f$ or \f$\nu_0=\mu_0,...,\nu_k>\mu_k\f$ for some k. Return 1 if true and 0 if equality and -1 otherwise. For example \f$(0,1,1,3,2)<_l(0,1,1,2,2)\f$ as abbcc would appear prior to abbdc in a dictionnary.
 
 * @param mi1 the first multi-index.
 * @param mi2 the second multi-index.
 * @return the boolean value
 */
int mpa_multi_index_lexical_less (const mpa_multi_index * mi1,const mpa_multi_index * mi2);

/*!
 Return the position of the active index that is <= i, result is equal to 1 if i is an active index and 0 otherwise.
 For exapmle if \f$\nu=(1,0,3,0,5,0)\f$, then the order of 2 is 1 and the order of 3 is also 1, but the result is 1 for the first and 0 for the second.
 *@param mi the multi-index.
 *@param i the index to search for.
 *@return a pair index-result.
 */
mpa_search_pair mpa_multi_index_order(const mpa_multi_index* mi,size_t i);

/*!
 * Creates the multi-index \f$\nu +e_i\f$ where \f$e_i=(0,0,...,1,0,...)\f$, returning a pointer to a newly initialized multi-index struct.
 * @param mi the multi-index .
 * @param i the index to increment.
 * @return pointer the new allocated multi-index.
 */
mpa_multi_index *mpa_multi_index_plus_ei(const mpa_multi_index* mi,size_t i);

/*!
* Creates the multi-index \f$\nu-e_i\f$ where \f$e_i=(0,0,...,1,0,...)\f$, returning a pointer to a newly initialized multi-index struct.
 An error is invoked if i is not active in the multi-index.
* @param mi the multi-index .
* @param i the index to decrement.
* @return pointer the new allocated multi-index.
*/
mpa_multi_index *mpa_multi_index_minus_ei(const mpa_multi_index* mi,size_t i);

#pragma mark -
#pragma mark functions to be preferabely used with a helper.

/*!
 Check if the \f$\mu=\nu+e_i\f$. Here \f$\mu\f$ and \f$\nu\f$ plays the role of mis and mi  The two multi-indices must have the same length
 and i less than length. Return 1 if true and 0 if not.
 * @param mi the first multi-index.
 * @param mis the second multi-index.
 * @param i the index
 * @return the boolean value
 */
int mpa_multi_index_isSuccessor(const mpa_multi_index* mi,const mpa_multi_index* mis,size_t i);

/*!
 Check if  \f$\mu=\nu-e_i\f$. Here \f$\mu\f$ and \f$\nu\f$ plays the role of mip and mi. The two multi-indices must have the same length and
 i less than length. Return 1 if true and 0 if not.
 * @param mi the first multi-index.
 * @param mip the second multi-index.
 * @param i the index
 * @return the boolean value
 */
int mpa_multi_index_isPredecessor(const mpa_multi_index* mi,const mpa_multi_index* mip,size_t i);

/*!
 Check if \f$\mu=\nu+e_i\f$. Here \f$\mu\f$ and \f$\nu\f$ plays the role of mis and mi. The two multi-indices must have the same length and
 i less than length. The parameter helper is an array of size length and all it's element must be equal to zero. Return 1 if true and 0 if not.
 * @param mi the first multi-index.
 * @param mis the second multi-index.
 * @param i the index
 * @param helper the array to use as helper.
 * @return the boolean value
 */
int mpa_multi_index_isSuccessor_helper(const mpa_multi_index* mi,const mpa_multi_index* mis,size_t i,size_t* helper);

/*!
 Check if \f$\mu=\nu-e_i\f$. her \f$\mu\f$ and \f$\nu\f$ plays the role of mip and mi. The two multi-indices must have the same length and
 i less than length. The parameter helper is an array of size length and all it's element must be equal to zero. Return 1 if true and 0 if not.
 * @param mi the first multi-index.
 * @param mip the second multi-index.
 * @param i the index
 * @param helper the array to use as helper.
 * @return the boolean value
 */
int mpa_multi_index_isPredecessor_helper(const mpa_multi_index* mi,const mpa_multi_index* mip,size_t i,size_t* helper);


#pragma mark -
#pragma mark function that where replaced by define, see mpa_multi_index.c

/*!
 Check if two index-coordinate pairs are element-wise equal, \f$(i_1,c_1)==(i_1,c_2)\f$ if and only if \f$i_1==i_2~and~c_1==c_2\f$.  Return 1 if true and 0 if not.
 * @param ic1 the first pair.
 * @param ic2 the second pair.
 * @return the boolean value
 */
	//int mpa_pairIC_equal(const mpa_pairIC * ic1,const mpa_pairIC * ic2);

/*!
 Compare two index-coordinate pairs in term of indices, \f$(i_1,c_1)<(i_1,c_2)\f$ if and only if \f$i_1<i_2\f$.
 Return >0 if true, 0 if equality and srictly smaller than 0 otherwise.
 * @param ic1 the first pair.
 * @param ic2 the second pair.
 * @return the comparaison result
 */
	//int mpa_pairIC_equal_index(const void* pair1,const void* pair2);

/*!
 Compare the two index-coordinate pairs .  Return 1 if they have the same index value and the coordinate value of the second is equal to the coordinate value of the second, 0 otherwise.
 * @param ic1 the first pair.
 * @param ic2 the second pair.
 * @return the boolean value.
 */
	//int mpa_pairIC_isSuccessor(const mpa_pairIC * ic1,const mpa_pairIC * ic2);

/*!
 Compare the two index-coordinates pair. \f$(i_1,c_1)<_l(i_1,c_2)\f$ if and only if \f$i_1<i_2\f$ or \f$i_1==i_2\f$ and \f$c_1>c_2\f$.
 Return 1 if ic1 lexically less than ic2, 0 if they are equal and -1 otherwise
 * @param ic1 the first pair.
 * @param ic2 the second pair.
 * @return the comparaison result.
 */
	//int mpa_pairIC_lexical_less(const mpa_pairIC * ic1,const mpa_pairIC * ic2);
/*!
 Return the factorial of a multi-index
*/


#pragma mark -
#pragma mark tensor product type function


/*!
 * Compute the factorial of a multi-index, which is the product of its coordinates factorials.
 * @param mi the multi-index .
 * @return the factorial
 */
size_t mpa_multi_index_factorial(const mpa_multi_index* mi);


/*!
 * Compute the power of vector by a multi-index, which is the product of coordinate-wise powers of
 * the vector by the multiindex.
 * @param x the vector.
 * @param mi the multi-index .
 * @return the power
 */
double mpa_multi_index_pow(const double* x, const mpa_multi_index* mi);


/*!
 * set the the term indexed by the multi-index of the Cartesian product sequence obtained from a given sequence. The
 * result is a newly allocated table of the same length as the multi-index, and set as prescribed
 * @param mi the multi-index .
 * @param seq the sequence r_0,r_1,r_2,.....
 * @return the sequence term associated
 */
double* mpa_multi_index_sequence_term(const double* seq, const mpa_multi_index* mi);


/*!
 * compute the entry indexed by a  multi-indices in the Cartesian product of a sequence.
 * @param array the sequence
 * @param nu the multi-index
 * @return the value of the product of array_{nu_j} for j=0,...,d-1
 */

double tensorProductEntry_1d (const double* array, const mpa_multi_index* nu);

/*!
 * compute the entry indexed by a pair of multi-indices in the Cartesian product of a matrix.
 * @param array the matrix
 * @param mu the first multi-index
 * @param nu the second multi-index
 * @return the value of the product of array_{mu_j, nu_j} for j=0,...,d-1
 */

double tensorProductEntry_2d (const double** array, const mpa_multi_index* mu, const mpa_multi_index* nu);


#endif /*__mpa_multi_index_h__*/


