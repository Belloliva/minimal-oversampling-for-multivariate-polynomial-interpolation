/*
 *  mpa_stdlib.h
 *  mpa
 *
 *  Created by Abdellah CHKIFA on 25/02/11.
 *  Copyright 2011 Oliva Electronics. All rights reserved.
 *
 */

#ifndef __mpa_stdlib_h__
#define __mpa_stdlib_h__


#include <stdlib.h>

/*! \struct mpa_search_pair
 \brief pair (index,result).
 
 A pair of two int used to get search result.
 */
typedef struct{
    int index;/*!< an int to hold the index of a search, it take value >=-1*/
    int result;/*!< boolean to hold the reslut of a search, it can be 0 or 1*/
}mpa_search_pair;


/*! \struct mpa_range
 \brief pair(position,length).
 
 A pair of two unsigned int used to apply functions to subarray.
 */
typedef struct
{
    size_t position;/*!< the start position of the range */
    size_t length;  /*!< the length of the range */
}
mpa_range;

/*!
 * Make a range from two positive integer.
 * @param p the range position.
 * @param l the range length.
 * @return the range of position p and length l
 */
mpa_range mpa_range_make(size_t p,size_t l);

/*!
 * Make a search pair from two integer.
 * @param ix the index of the search pair.
 * @param rt the result of the search pair, 1 or 0.
 * @return the search pair of index p and result rt l
 */
mpa_search_pair sp_make(int ix,int rt);

/*!
 * Swap the content of two pointer.
 * @param a the first pointer.
 * @param b the second pointer.
 * @param size the size in bytes of the what is being swapped.
 * @return EXIT_SUCCESS is everything went well and EXIT_FAILURE if not
 */
int mpa_swap (void *const a, void *const b, size_t size);

/*!
 * Perform a binary search in an array and find the index of.
 * @param key ....
 * @param base .
 * @param nel  the number of the element of base
 * @param width the size in bytes of the element in the array
 * @param compar the function used in comparaion
 * @return {ix,1} if key was found in the array and ix the index corresponding to key or {ix,0} if key was not found,
 * ix is the index of element immeditely before key, so if we insert key in base, its index will be ix+1
 */
mpa_search_pair mpa_bsearch(const void *key, const void *base, size_t nel,
                            size_t width, int (*compar)(const void *, const void *));

void mpa_insert(const void *key, const void *baseSrc, void *baseDest, size_t nel, size_t index, size_t width);
void mpa_delete                 (const void *baseSrc, void *baseDest, size_t nel, size_t index, size_t width);

int mpa_diff(const void *base1, size_t nel1,
             const void *base2, size_t nel2,
             size_t width, int (*compar)(const void *, const void *));

void mpa_delete_append(void *base, size_t nel,size_t index,size_t width);

#endif /*__mpa_stdlib_h__*/
