/*
 *  mpa_multi_index.c
 *  bulk_chase
 *
 *  Created by Abdellah CHKIFA on 17/02/11.
 *  Copyright 2011 oliva Electronics. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "mpa_defines.h"
#include "mpa_multi_index.h"
#include "mpa_combinatorics.h"

/*!
 Check if two index-coordinate pairs are element-wise equal, \f$(i_1,c_1)==(i_1,c_2)\f$ if and only if \f$i_1==i_2~and~c_1==c_2\f$.  Return 1 if true and 0 if not.
 * @param ic1 the first pair.
 * @param ic2 the second pair.
 * @return the boolean value
 */

#define mpa_pairIC_equal(ic1,ic2)\
((ic1)!=NULL && (ic2)!=NULL && (ic1)->index==(ic2)->index && (ic1)->coordinate == (ic2)->coordinate)

/*!
 Compare the two index-coordinate pairs, return 1 if they have the same index value and the coordinate value of the second is equal to the coordinate
 value of the first plus one, 0 otherwise.
 * @param ic1 the first pair.
 * @param ic2 the second pair.
 * @return the boolean value.
 */

#define mpa_pairIC_successor(ic1,ic2)\
((ic1)!=NULL && (ic2)!=NULL && (ic1)->index==(ic2)->index && ((ic1)->coordinate+1) == (ic2)->coordinate)

/*!
 Compare the two index-coordinate pairs, i.e. \f$(i_1,c_1)<_l(i_1,c_2)\f$ if and only if \f$i_1<i_2\f$ or \f$i_1==i_2\f$ and \f$c_1>=c_2\f$,
 return >0 if ic1 lexically less than ic2, 0 if they are equal and <0 otherwise
 * @param ic1 the first pair.
 * @param ic2 the second pair.
 * @return the comparaison result.
 */
/*
#define mpa_pairIC_lexical_less(ic1,ic2)\
((ic2)->index - (ic1)->index) + ((ic2)->index == (ic1)->index) * ((ic1)->coordinate - (ic2)->coordinate)
*/

static int mpa_pairIC_lexical_less(const void* pair1,const void* pair2)
{
    mpa_pairIC* ic1 = (mpa_pairIC*)pair1;
    mpa_pairIC* ic2 = (mpa_pairIC*)pair2;
    int equl= (int)ic2->index - (int)ic1->index;
    if (equl !=0)
        return equl;
    else
        return
        (int)ic1->coordinate - (int)ic2->coordinate;
}


/*!
 Compare two index-coordinate pairs in term of indices, \f$(i_1,c_1)<(i_1,c_2)\f$ if and only if \f$i_1<i_2\f$,
 return >0 if true, 0 if equality and <0 otherwise.
 * @param ic1 the first pair.
 * @param ic2 the second pair.
 * @return the comparaison result
 */

//#define mpa_pairIC_index_less(ic1,ic2) ((ic2)->index - (ic1)->index)


static int mpa_pairIC_index_less(const void* pair1,const void* pair2)
{return (int)(((mpa_pairIC*)pair1)->index - ((mpa_pairIC*)pair2)->index);}


#pragma mark -
#pragma mark memroy management

mpa_multi_index *mpa_multi_index_alloc (const size_t n,const size_t l)
{
    if(n>l)
        perror("the multi_index size must be less than the multi_index length\n");
    
    mpa_multi_index* mi =(mpa_multi_index*) malloc(sizeof(mpa_multi_index));
    if (mi == NULL){
        perror("failed to allocate memory for multi-index struct.\n");
        return NULL;
    }
    else{
        if (n==0) {//the null multi-index has all its entries zeros so has no data
            mi->size=n;
            mi->length=l;
            mi->data=NULL;
        }
        else{
            mpa_pairIC* ic =(mpa_pairIC*) calloc(n,sizeof(mpa_pairIC));
            if (ic == NULL)
                perror("failed to allocate memory for multi-index datas.\n");
            else {
                mi->size=n;
                mi->length=l;
                mi->data=ic;
            }
        }
        return mi;
    }
}

void mpa_multi_index_memcpy(mpa_multi_index * dest,const mpa_multi_index * src)
{
    if(dest!=NULL && src!=NULL && dest->size != src->size) // || dest->length!=src->length)
        perror("cannot copy to a multi-index that has a different size.\n");
    mpa_pairIC *icd = dest->data, *ics = src->data;
    size_t j=src->size;
    while (j--)
        *icd++=*ics++;
}

void mpa_multi_index_free(mpa_multi_index* mi)
{
    free(mi->data);
    free(mi);
}

#pragma mark -
#pragma mark getters

mpa_pairIC* mpa_multi_index_get  (const mpa_multi_index * mi,size_t i)
{
    if (i>=mi->length)
        perror("cannot get the index-coordinate pair, the index is out of range");
    mpa_search_pair osp = mpa_multi_index_order(mi, i);
    return (osp.result)?&mi->data[osp.index]:NULL;
}

size_t* mpa_multi_index_vector(const mpa_multi_index * mi)
{
    size_t* coorV = (size_t*) calloc(mi->length, sizeof(size_t));
    if (coorV == NULL)
        perror("failed to allocate memory for multi-index vector representation.");
    else{
        for (int i=0; i<mi->length; i++) { coorV[i]=0;}
        mpa_pairIC *ic=mi->data;
        for (int i=0; i<mi->size; i++) {
            coorV[ic[i].index]=ic[i].coordinate;
        }
    }
    return coorV;
}

#pragma mark -
#pragma mark simple functions

/*!
 The l1 norm of the multi-index.
 * @param mi the multiindex.
 * @return the l1 norm
 */
size_t mpa_multi_index_abs(const mpa_multi_index * mi)
{
    size_t absV=0,j=mi->size;
    mpa_pairIC* it=mi->data;
    while (j--)
        absV+=it++->coordinate;
    return absV;
}

/*!
 The inner product of a multi-index with a vector of weights.
 * @param mi the multiindex.
 * @param weight the vector of wiehgts.
 * @return the inner product
 */
double mpa_multi_index_dot_vector(const mpa_multi_index * mi,const double* weight)
{
    double sum=0.0;
    size_t	j=mi->size;
    mpa_pairIC* it=mi->data;
    while (j--){
        sum+=it->coordinate * weight[it->index];
        it++;
    }
    return sum;
}


/*!
 The maximum value of the multi-index coordinates .
 * @param mi the multiindex.
 * @return the maximum coordinate
 */
size_t mpa_multi_index_max(const mpa_multi_index * mi)
{
    size_t max=0,j=mi->size,maxV;
    mpa_pairIC* it=mi->data;
    while (j--)
        if ((maxV=it++->coordinate)>max)
            max=maxV;
    return max;
    /*you should not be tempted by a ternary expression because
     it replaces max by max when it->coordinate<max which is useless.*/
}

/*!
 The minimum value of the multi-index coordinates .
 * @param mi the multiindex.
 * @return the minimum coordinate
 */
size_t mpa_multi_index_min(const mpa_multi_index * mi)
{
    if (mi->size<mi->length)
        return 0;
    mpa_pairIC* it=mi->data;
    size_t min=it++->coordinate,j=mi->size-1,minV;
    while (j--)
        if ((minV=it++->coordinate)<min)
            min=minV;
    return min;
}

#pragma mark -
#pragma mark i/o functions

/*!
 Print all the coordinates of the multi-index as a vector.
 * @param mi the multiindex.
 */
void mpa_multi_index_printf(const mpa_multi_index * mi)
{
    size_t* data=mpa_multi_index_vector(mi);
    printf("(");
    size_t i;
    for (i=0; i<mi->length-1; i++)
        printf("%lu,",data[i]);
    printf("%lu)\n",data[i]);
    free(data);
}


/*!
 Print all the coordinates of the multi-index to a file.
 * @param fp the file.
 * @param mi the multiindex.
 */

void mpa_multi_index_fprintf(FILE *fp,const mpa_multi_index* mi)
{
    size_t i=0,j=mi->size;
    mpa_pairIC* it=mi->data;
    char buffer[5*mi->length];
    char* buff=buffer;
    
    while (i<mi->length) {
        if (j && it->index==i ) {
            buff+=sprintf(buff,"%lu ",it++->coordinate);
            j--;
        }
        else
            buff+=sprintf(buff,"%d ",0);
        i++;
    }
    fprintf(fp,"%s \n",buffer);
}


#pragma mark -
#pragma mark Comparaison functions


/*!
 check if two multi-indices are element-wise equal.  Return 1 if true and 0 if not.
 * @param mi1 the first multi-index.
 * @param mi2 the second multi-index.
 * @return the boolean value
 */

int mpa_multi_index_equal(const mpa_multi_index * mi1,const mpa_multi_index * mi2)
{
    if (mi1->size!=mi2->size)// || mi1->length!=mi2->length)
        return 0;
    mpa_pairIC *it1=mi1->data,*it2=mi2->data;
    size_t j=0;
    while (j<mi1->size && mpa_pairIC_equal(it1++, it2++))
        j++;
    return (j==mi1->size);
}


/*!
 check if two multi-indices are element-wise less than or equal.  Return 1 if true and 0 if not.
 * @param mi the first multi-index.
 * @param vec the second multi-index but in vector format.
 * @return the boolean value
 */

int mpa_multi_index_less_vector (const mpa_multi_index * mi,const size_t* vec)
{
    for(int j=0;j<mi->size;j++){
        if (mi->data[j].coordinate > vec[mi->data[j].index])
            return 0;
    }
    return 1;
}


/*!
 check if a multi-index is element-wise less than and is not equal to another multi-index.  Return 1 if true and 0 if not.
 * @param mi1 the first multi-index.
 * @param mi2 the second multi-index.
 * @return the boolean value
 */

int mpa_multi_index_lessOrEqual (const mpa_multi_index * mi1,const mpa_multi_index * mi2)
{
    if (mi1->size>mi2->size)
        return 0;
    
    size_t mi2_vec[mi2->length];
    for (int k=0; k<mi2->length; k++)
        mi2_vec[k]=0;
    
    for (int k=0; k<mi2->size ; k++)
        mi2_vec[mi2->data[k].index]=mi2->data[k].coordinate;
    
    for (int k=0; k<mi1->size ; k++) {
        size_t ix_j = mi1->data[k].index;
        size_t nu_j = mi1->data[k].coordinate;
        if (mi2_vec[ix_j] < nu_j)
            return 0;
    }
    return 1;
    /*
    int diff  =0;
    int count =0;
    mpa_pairIC* ic2 = mi2->data;
    for(int j=0;j<mi1->size;j++){
        size_t ix_j = mi1->data[j].index;
        size_t nu_j = mi1->data[j].coordinate;
        while (ic2->index < ix_j) {ic2++;count++;}
        if (ic2->index == ix_j && nu_j <= ic2->coordinate)
            diff += ic2->coordinate - nu_j;
        else
            return 0;
    }
    return diff >0 || count < mi2->size;
     */
}


int mpa_multi_index_less (const mpa_multi_index * mi1,const mpa_multi_index * mi2)
{
    if (mi1->size>mi2->size)
        return 0;
    
    size_t mi1_vec[mi1->length];
    size_t mi2_vec[mi2->length];
    for (int k=0; k<mi1->length; k++){
        mi1_vec[k]=0;
        mi2_vec[k]=0;
    }
    
    for (int k=0; k<mi1->size ; k++){
        size_t i    = mi1->data[k].index;
        size_t nu_i = mi1->data[k].coordinate;
        mi1_vec[i]=nu_i;
    }
    
    for (int k=0; k<mi2->size ; k++){
        size_t i    = mi2->data[k].index;
        size_t mu_i = mi2->data[k].coordinate;
        mi2_vec[i]=mu_i;
    }
    
    int diff = 0;
    
    for (int k=0; k<mi1->length ; k++) {
        size_t nu_k = mi1_vec[k];
        size_t mu_k = mi2_vec[k];
        if (mu_k < nu_k)
            return 0;
        diff+= (mu_k - nu_k);
    }
    /*
    for (int k=0; k<mi1->size ; k++) {
        size_t ix_j = mi1->data[k].index;
        size_t nu_j = mi1->data[k].coordinate;
        if (mi2_vec[ix_j] < nu_j)
            return 0;
        diff+= mi2_vec[ix_j] - nu_j;
    }*/
    return (diff>0);
    /*
    int diff  =0;
    int count =0;
    mpa_pairIC* ic2 = mi2->data;
    for(int j=0;j<mi1->size;j++){
        size_t ix_j = mi1->data[j].index;
        size_t nu_j = mi1->data[j].coordinate;
        while (ic2->index < ix_j) {ic2++;count++;}
        if (ic2->index == ix_j && nu_j <= ic2->coordinate)
            diff += ic2->coordinate - nu_j;
        else
            return 0;
    }
    return diff >0 || count < mi2->size;
     */
}


/*!
 check if two multi-indices are element-wise less than or equal.  Return 1 if true and 0 if not.
 * @param mi1 the first multi-index.
 * @param mi2 the second multi-index.
 * @param helper an auxiliary table of the same length used in the comparaison process.
 * @return the boolean value
 */

int mpa_multi_index_less_helper (const mpa_multi_index * mi1,
                                 const mpa_multi_index * mi2,
                                 size_t* helper)
{
    if (mi1->size>mi2->size)
        return 0;
    
    for (int k=0; k<mi1->length; k++) {helper[k]=0;}
    for (int k=0; k<mi1->size ; k++) {
        helper[mi1->data[k].index]=mi1->data[k].coordinate;
    }
    /*
    mpa_pairIC *ic1=mi1->data;
    size_t j=mi1->size;
    while (j--) {
        helper[ic1->index]=ic1->coordinate;
        ic1++;
    }
    */
    
    int result=0,diff;
    mpa_pairIC *ic2=mi2->data;
    int j=0;
    while (j<mi2->size && (diff=(int)ic2->coordinate - (int)helper[ic2->index])>=0) {
        result|=(diff>0);
        ic2++;
        j++;
    }
    return (j==mi2->size && result);
}


/*!
 Compare two multi-indices using lexical compaison.  Return 1 if mi1 is lexically less than mi2, 0 if they are equal and -1 otherwise
 * @param mi1 the first multi-index.
 * @param mi2 the second multi-index.
 * @return the boolean value
 */

int mpa_multi_index_lexical_less (const mpa_multi_index * mi1,const mpa_multi_index * mi2)
{
    mpa_pairIC *ic1=mi1->data,*ic2=mi2->data;
    size_t j=0, minSize=(mi1->size<mi2->size)?mi1->size:mi2->size;
    while (j<minSize && mpa_pairIC_equal(ic1,ic2)){
        j++;
        ic1++;
        ic2++;
    }
    if (j==minSize)
        return (int)(mi2->size-mi1->size);
    return mpa_pairIC_lexical_less(ic1, ic2);
}

/*!
 find the pair index-coordinate corresponding to an integer i in a multiiindex.
 * @param mi the multi-index.
 * @param i the integer .
 */


mpa_search_pair mpa_multi_index_order(const mpa_multi_index* mi,size_t i)
{
    if (i>=mi->length)
        perror("index is out of range, index must be less than the multi-index length");
    
    if (mi->data==NULL)
        return sp_make(-1,0);
    
    mpa_pairIC toSearch4={i,0};//0 here is unused
    return mpa_bsearch(&toSearch4, mi->data, mi->size, sizeof(mpa_pairIC), mpa_pairIC_index_less);
}


#pragma mark -
#pragma mark successors and predecessors


/*!
 construct the multi-index which results from adding a kronecker vector ei to a multiiindex.
 * @param mi the multi-index.
 * @param i the index of the kronecker vector.
 */

mpa_multi_index *mpa_multi_index_plus_ei(const mpa_multi_index* mi,size_t i)
{
    mpa_search_pair order=mpa_multi_index_order(mi,i);
    if (order.result==0){
        mpa_multi_index* mis=mpa_multi_index_alloc(mi->size+1,mi->length);
        mpa_pairIC pair = {i,1};
        mpa_insert(&pair, mi->data, mis->data, mi->size,order.index+1, sizeof(mpa_pairIC));
        return mis;
    }
    else{
        mpa_multi_index* mis=mpa_multi_index_alloc(mi->size,mi->length);
        mpa_multi_index_memcpy(mis,mi);
        mis->data[order.index].coordinate+=1;
        return mis;
    }
}

/*!
 construct the multi-index which result from substracting a kronecker vector ei to a multiiindex.
 * @param mi the multi-index.
 * @param i the index of the kronecker vector.
 */

mpa_multi_index *mpa_multi_index_minus_ei(const mpa_multi_index* mi,size_t i)
{
    
    mpa_search_pair order = mpa_multi_index_order(mi,i);
    if (order.result==0){
        perror("cannot decrement an inactive coordinate");
        return NULL;
    }
    else{
        if (mi->data[order.index].coordinate>1){
            mpa_multi_index* mip=mpa_multi_index_alloc(mi->size,mi->length);
            mpa_multi_index_memcpy(mip,mi);
            mip->data[order.index].coordinate--;
            return mip;
        }
        else{
            mpa_multi_index* mip=mpa_multi_index_alloc(mi->size-1,mi->length);
            if (mi->size==1)
                return mip;//mip=0 the null multi-index in this case
            mpa_delete(mi->data, mip->data, mi->size, order.index, sizeof(mpa_pairIC));
            return mip;
        }
    }
}

/*!
 say if a multi-index is the successor of a multi-index, meaning it is equal to it plus a kronecker vector. Return 1 if true and 0 if not
 * @param mi the first multi-index.
 * @param mis the second multi-index.
 * @param i the index of the kronecker vector.
 * @return the boolean value
 */


int mpa_multi_index_isSuccessor(const mpa_multi_index* mi,const mpa_multi_index* mis,size_t i)
{
    //	if (mi->length != mis->length)
    //perror("the two multi-indices must have the same length");
    int diff = (int) (mis->size-mi->size);
    if (diff<0 || diff>1)
        return 0;
    mpa_pairIC* it = mi->data;
    mpa_pairIC* its=mis->data;
    int j=0;
    int iExist=0;
    
    if(diff==0){
        while (j++<mi->size){
            if( it->index!=i && !mpa_pairIC_equal(it,its))
                return 0;
            else if( its->index==i && (iExist=1) && !mpa_pairIC_equal(it,its))// iExist assignemnt not check!!!
                return 0;
            it++;
            its++;
        }
    }
    
    else if(diff==1){
        while (j++<mi->size){
            if(its->index!=i && !mpa_pairIC_equal(it,its))
                return 0;
            else if(its->index==i && (iExist=1)){// iExist assignemnt not check!!!
                if(its->coordinate!=1)
                    return 0;
                its++;//we advance to compare equalities again for indices >i
            }
            else {
                it++;
                its++;
            }
        }
        if (iExist==0)//we have exhausted mi and we found equalities so far
            iExist = (its->index==i && its->coordinate==1);
    }
    return iExist;
}

/*!
 say if a multi-index is the predecessor of a multi-index, meaning it is equal to it minus a kronecker vector. Return 1 if true and 0 if not
 * @param mi the first multi-index.
 * @param mip the second multi-index.
 * @param i the index of the kronecker vector.
 * @return the boolean value
 */


int mpa_multi_index_isPredecessor(const mpa_multi_index* mi,const mpa_multi_index* mip,size_t i)
{
    return mpa_multi_index_isSuccessor(mip,mi ,i);
}


/*!
 say if a multi-index is the successor of a multi-index, meaning it is equal to it plus a kronecker vector, with the use of a helper. Return 1 if true and 0 if not
 * @param mi the first multi-index.
 * @param mis the second multi-index.
 * @param i the index of the kronecker vector.
 * @param helper the array used as helper.
 * @return the boolean value
 */

int mpa_multi_index_isSuccessor_helper(const mpa_multi_index* mi,const mpa_multi_index* mis,size_t i,size_t* helper)
{
    //	if (mi->length != mis->length)
    //	perror("the two multi-indices must have the same length");
    int diff = (int)(mis->size-mi->size);
    if (diff<0 || diff>1)
        return 0;
    
    for (int j=0; j<mi->size; j++) {
        helper[mi->data[j].index]=mi->data[j].coordinate;
    }
    
    mpa_pairIC* it=mis->data;
    
    for (int j=0; j<mis->size; j++) {
        if (helper[it->index]==0)
            return 0;
        else {
            if (it->index==i && it->coordinate!=(helper[i]+1))
                return 0;
            else if(it->coordinate!=helper[it->index])
                return 0;
        }
        it++;
    }
    return 1;
    
}


/*!
 say if a multi-index is the predessecor of a multi-index, meaning it is equal to it minus a kronecker vector, with the use of a helper. Return 1 if true and 0 if not
 * @param mi the first multi-index.
 * @param mip the second multi-index.
 * @param i the index of the kronecker vector.
 * @param helper the array used as helper.
 * @return the boolean value
 */

int mpa_multi_index_isPredecessor_helper(const mpa_multi_index* mi,const mpa_multi_index* mip,size_t i,size_t* helper)
{
    return mpa_multi_index_isSuccessor_helper(mip,mi ,i ,helper);
}

#pragma mark -
#pragma mark functions that give boolean values


int mpa_multi_index_isZero(const mpa_multi_index * mi)
{
    if (mi->size ==0)
        return True;
    else
        return False;
}


int mpa_multi_index_isUnitVector(const mpa_multi_index * mi)
{
    if (mi->size ==1 && mi->data[0].coordinate ==1)
        return True;
    else
        return False;
}

/*!
 say if a multi-index is within a hypercube of a give side. Return 1 if true and 0 if not
 * @param mi the first multi-index.
 * @param side the side of the hypercube.
 * @return the boolean value
 */
/* this is equivalent to max(mi)<= side but max require to iterate throught
 mi while here we could get false by only checking the first coordinate*/
int mpa_multi_index_isInsideHypercube(const mpa_multi_index * mi,size_t side)
{
    mpa_pairIC* it = mi->data;
    size_t j=mi->size;
    while (j--)
        if (it++->coordinate>side)
            return 0;
    return 1;
}

/*!
 say if a multi-index is within a simplex of a give shape. Return 1 if true and 0 if not
 * @param mi the first multi-index.
 * @param weight the weights used in the definition of the simplex.
 * @param k the width used in the definition of the simplex.
 * @return the boolean value
 */


int mpa_multi_index_isInsideSimplex(const mpa_multi_index * mi,double* weight,double k)
{
    return (mpa_multi_index_dot_vector(mi,weight)<=k);
}

#pragma mark -
#pragma mark functions as products

/*!
 compute the factorial of a multi-index.
 * @param mi the multi-index.
 * @return the factorial value
 */

size_t mpa_multi_index_factorial(const mpa_multi_index* mi)
{
    size_t fact=1,i=mi->size;
    mpa_pairIC *ic=mi->data;
    while (i--)
        fact*=mpa_factorial(ic++->coordinate);
    return fact;
}


/*!
 compute the power of a vector by a multi-index.
 * @param x the vector .
 * @param mi the multi-index.
 * @return the power value
 */
double mpa_multi_index_pow(const double* x, const mpa_multi_index* mi)
{
    double prod=1.;
    size_t j=mi->size;
    mpa_pairIC *ic=mi->data;
    while (j--){
        prod*=pow(x[ic->index],ic->coordinate);
        ic++;
    }
    return prod;
}

double* mpa_multi_index_sequence_term(const double* seq, const mpa_multi_index* mi)
{
    double* seq_mi = (double*) calloc(mi->length, sizeof(double));
    
    for (size_t j=0; j<mi->length; j++)
        seq_mi[j] = seq[0];
    
    for (size_t j=0; j<mi->size; j++){
        size_t i    = mi->data[j].index;
        size_t nu_i = mi->data[j].coordinate;
        seq_mi[i] = seq[nu_i];
    }
    return seq_mi;
}


double tensorProductEntry_1d (const double* array, const mpa_multi_index* nu)
{
    double result = 1.;
    for (int j=0; j<nu->size; j++) {
        //size_t i    = nu->data[j].index;
        size_t nu_i = nu->data[j].coordinate;
        result     *= array[nu_i];
    }
    return result;
}


double tensorProductEntry_2d (const double** array, const mpa_multi_index* nu, const mpa_multi_index* mu)
{
    double result = 1.;
    size_t* muVec = mpa_multi_index_vector(mu);
    for (int j=0; j<nu->size; j++) {
        size_t i    = nu->data[j].index;
        size_t nu_i = nu->data[j].coordinate;
        result     *= array[nu_i] [muVec[i]];
    }
    free(muVec);
    return result;
}
