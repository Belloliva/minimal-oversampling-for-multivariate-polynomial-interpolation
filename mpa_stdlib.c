/*
 *  mpa_stdlib.c
 *  mpa
 *
 *  Created by Abdellah CHKIFA on 25/02/11.
 *  Copyright 2011 Oliva Electronics. All rights reserved.
 *
 */
#include "mpa_stdlib.h"
#include <string.h>


int mpa_swap (void *const a, void *const b, size_t size)
{
    void *const temp = malloc(size);
    if ( temp == NULL )
        return EXIT_FAILURE;
    
    memcpy(temp, a, size);
    memcpy(a, b, size);
    memcpy(b, temp, size);
    
    free(temp);
    return EXIT_SUCCESS;
}

// simply make a range
mpa_range mpa_range_make(size_t p,size_t l)
{
    mpa_range R={p,l};
    return R;
}

// simply make a search pair
mpa_search_pair sp_make(int in, int rt)
{
    mpa_search_pair sp={in,rt};
    return sp;
}
 
mpa_search_pair mpa_bsearch(const void *key,
                            const void *base,
                            size_t nel,
                            size_t width,
                            int (*compar)(const void *, const void *))
{
    int left=0,right=(int)nel-1,cursor;
    int cmpl,cmpr;//compaison left and right.
    
    cmpr = compar(key,base+right*width);
    if (cmpr>=0)
        return (cmpr==0)?sp_make(right,1):sp_make(right,0);
    
    cmpl = compar(key,base+left*width);
    if (cmpl<=0)
        return (cmpl==0)?sp_make(0,1):sp_make(-1,0);
    
    int cmp;
    while (right-left>1) {
        cursor=(left+right)/2;
        cmp= compar(key,base+cursor*width);
        if (cmp>0)
            left=cursor;
        else if (cmp<0)
            right=cursor;
        else
            return sp_make(cursor, 1);
    }
    return sp_make(left,0);
}

void mpa_insert(const void *key, const void *baseSrc,void *baseDest, size_t nel, size_t index, size_t width)
{
    if (index==0){
        memmove(baseDest+width,baseSrc,width*nel);
        memcpy(baseDest,key, width);
    }
    /*make sure every time it is not the same array, otherwise, we perform a self
     copy that may slow the code for large arrays*/
    else if (index==nel) {
        if (baseSrc!=baseDest)
            memmove(baseDest, baseSrc, width*nel);
        memcpy(baseDest+(nel*width), key,width);
    }
    
    else{
        if (baseSrc!=baseDest)
            memmove(baseDest, baseSrc, width*index);
        memmove(baseDest+(index+1)*width, baseSrc+index*width, width*(nel-index));
        memcpy(baseDest+index*width,key,width);//dont mess up the order
    }
}

void mpa_delete(const void *baseSrc,void *baseDest, size_t nel,size_t index,
                size_t width)
{
    if (index==0){
        for (int i=0; i<nel-1;i++)
            memcpy(baseDest+i*width, baseSrc+(i+1)*width, width);
    }
    else if (index==nel-1) {
        //make sure it's not the same array, otherwise, we perform a self
        //copy that may slow the code for large arrays
        if (baseSrc!=baseDest) {
            for (int i=0; i<nel-1;i++)
                memcpy(baseDest+i*width, baseSrc+i*width, width);
        }
    }
    else{
        for (int i=0; i<index;i++)
            memcpy(baseDest+i*width, baseSrc+i*width, width);
        for (int i=(int)index; i<nel-1;i++)
            memcpy(baseDest+i*width, baseSrc+(i+1)*width, width);
    }
}

int mpa_diff(const void *base1, size_t nel1,
             const void *base2, size_t nel2,
             size_t width, int (*compar)(const void *, const void *))
{
    if (nel1==0 || nel2==0)
        return 0;
    if(nel2<nel1)
        return mpa_diff(base2,nel2,base1,nel1,width,compar);
    int cmp;
    if ((cmp=compar(base1+(nel1-1)*width,base2))<=0)
        return (cmp==0)?1:0;
    if ((cmp=compar(base1,base2+(nel2-1)*width))>=0)
        return (cmp==0)?1:0;
    
    size_t size2=nel2,size1=nel1;
    int sum=0;
    mpa_search_pair osp;
    void *key=(void*)base1;
    void *begin=(void*)base2;
    
    while (size2>0&&size1>0) {
        osp=mpa_bsearch(key,begin,size2, width,compar);
        sum+=osp.result;
        if (osp.index!=-1){
            size2-=osp.index;
            begin+=osp.index*width;
        }
        key+=width;
        size1--;
    }
    return sum;
}


