/*
 *  mpa_monotone_graph.c
 *  bulk_chase
 *
 *  Created by Abdellah CHKIFA on 24/02/11.
 *  Copyright 2011 oliva Electronics. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>

#include "mpa_defines.h"
#include "mpa_monotone_graph.h"
#include "mpa_combinatorics.h"


mg_node* mg_node_alloc(mpa_multi_index* miPtr)
{
    mg_node* node = (mg_node*)malloc(sizeof(mg_node));
    node->IN=0;
    node->info  = NULL;
    node->miPtr =miPtr;
    node->mi_abs = mpa_multi_index_abs(miPtr);
    node->successors = (mg_node**)calloc(miPtr->length, sizeof(mg_node*));
    node->predecessors = (mg_node**)calloc(miPtr->size  , sizeof(mg_node*));
    return node;
}

void mg_node_free(mg_node* node)
{
    if (node) {
        /* every node is only responsible for freeing her multi-index*/
        mpa_multi_index_free(node->miPtr);
        // free(node->info); info is freed elsewhere because it is no more than a cast
        // to a void pointer and hence hold no information on how memory was allocated
        free(node->successors);//this only frees the successors array not the successors themselves
        free(node->predecessors);//this only frees the predecessors array not the predecessors themselves
        free(node);
    }
}

#pragma mark -
#pragma mark Compaison functions


int mg_node_compare_abs(const void* node1,const void* node2)
{
    //int c1= (int)mpa_multi_index_abs((*(mg_node**)node1)->miPtr);
    //int c2= (int)mpa_multi_index_abs((*(mg_node**)node2)->miPtr);
    //return c1-c2;

    int c1= (int)((*(mg_node**)node1)->mi_abs);
    int c2= (int)((*(mg_node**)node2)->mi_abs);
    return c1-c2;
}

#pragma mark -
#pragma mark setting functions

void mg_node_set_successor(mg_node* node, mg_node* succ,int i)
{
    node->successors[i]=succ;
    mpa_search_pair p = mpa_multi_index_order(succ->miPtr,i);
    succ->predecessors[p.index]=node;
    /*in the last line, don't call the function mg_node_set_predecessor
     to set the predecessor, because it will set succ again as successor to node */
}

void mg_node_set_predecessor(mg_node* node, mg_node* pred,int p)
{
    node->predecessors[p]=pred;
    size_t i = node->miPtr->data[p].index;
    pred->successors[i]=node;
}

void mg_node_set_all_predecessors(mg_node* node, mg_node* aPred,size_t p)
{
    if (node->miPtr->size==1){
        node->predecessors[0]=aPred;
        aPred->successors[p] =node;
    }
    else {
        size_t i,ip;
        // i will iterate on node predecessors and ip on aPred predecessors
        // we know that if i is in 0...n then ip is in 0...n or in 0...n-1 whether
        // the node multi-index coordinate at p is >1 or is 1.
        // this will be decided inside this loop.
        size_t cI; //current index
        mg_node* cNode; //current node
        for (i=0,ip=0; i<node->miPtr->size;i++,ip++) {
            if ((cI=node->miPtr->data[i].index)==p){
                node->predecessors[i]=aPred;
                cNode=aPred;
                if (node->miPtr->data[i].coordinate==1)
                    ip--;
            }
            else{
                cNode = aPred->predecessors[ip]->successors[p];
                node->predecessors[i]=cNode;
            }
            cNode->successors[cI]=node;
        }
    }
}

int  mg_node_all_predecessors_in(const mg_node* node)
{
    mpa_multi_index* mi = node->miPtr;
    int is_admissible =1;
    for (int p=0; p<mi->size; p++){
        mg_node* pred = node->predecessors[p];
        is_admissible &=(pred!=NULL && pred->IN==1);
    }
    return is_admissible;
}

int  mg_node_successor_all_predecessors_in(const mg_node* node, size_t k)
{
    mpa_multi_index* mi = node->miPtr;
    int is_admissible =1;
    for (int p=0; p<mi->size; p++){
        mg_node* pred_node = node->predecessors[p];
        mg_node* pred_succ = pred_node->successors[k];
        is_admissible &=(pred_succ!=NULL && pred_succ->IN==1);
    }
    return is_admissible;
}

mpa_monotone_graph* mpa_monotone_graph_alloc(size_t capacity)
{
    mpa_monotone_graph* omg=(mpa_monotone_graph*) malloc(sizeof(mpa_monotone_graph));
    omg->capacity = capacity;
    omg->size = 0;
    omg->nodes =(mg_node**) calloc(capacity, sizeof(mg_node*));
    return omg;
}

void mpa_monotone_graph_free(mpa_monotone_graph* omg)
{
    for(int i=0;i<omg->size;i++)
        mg_node_free(omg->nodes[i]);
    free(omg->nodes);
    free(omg);
}

void mpa_monotone_graph_resize(mpa_monotone_graph* omg,size_t newCapacity)
{
    if (newCapacity<omg->capacity){
        for (size_t i=newCapacity+1;i<omg->capacity;i++)
            mg_node_free(omg->nodes[i]);
    }
    omg->capacity=newCapacity;
    if (omg->size>newCapacity)
        omg->size=newCapacity;
    
    mg_node** new_nodes =(mg_node**) calloc(newCapacity, sizeof(mg_node*));
    for (int i=0; i<omg->size; i++)
        new_nodes[i]=omg->nodes[i];
    
    free(omg->nodes);
    omg->nodes=new_nodes;
}

void mpa_monotone_graph_add_node(mpa_monotone_graph* omg,mg_node* node)
{
    if (omg->size==omg->capacity)
        mpa_monotone_graph_resize(omg,2*omg->capacity);
    omg->nodes[omg->size]=node;
    omg->size++;
}

void mpa_monotone_graph_insert_node(mpa_monotone_graph* omg,
                                    mpa_range r,
                                    mg_node* node,
                                    int (*cmp)(const void *, const void *))
{
    if (r.position+r.length>=omg->capacity){
        mpa_monotone_graph_resize(omg,2*omg->capacity);
    }
    if (r.length==0){
        omg->nodes[r.position]=node;
        omg->size++;
        return;
    }
    mpa_search_pair osp = mpa_bsearch(&node, omg->nodes+r.position, r.length, sizeof(mg_node*), cmp);
    //mpa_insert_new(&node,omg->nodes+r.position,r.length, osp.index+1,sizeof(mg_node*));
    //mpa_insert_to_graph(node,omg->nodes+r.position, r.length,osp.index+1);
    mpa_insert(&node, omg->nodes+r.position, omg->nodes+r.position, r.length, osp.index+1,sizeof(mg_node*));
    omg->size++;
}



int mg_node_compare_suite(const void* node1,const void* node2)
{
    double res= *(double*)((*(mg_node**)node2)->info)- *(double*)((*(mg_node**)node1)->info);
    return signC(res);
}

mpa_monotone_graph* mpa_monotone_graph_envellope(size_t dim, size_t capacity,mpa_function func)
{
    double* d;
    mpa_multi_index* mi;
    mpa_multi_index* succ_mi;
    mg_node* node;
    mg_node* succ_node;

    mpa_monotone_graph* omg = mpa_monotone_graph_alloc(capacity);
    
    /*####### create and add the root node to the monotone graph #######*/
    mi = mpa_multi_index_alloc(0, dim);
    node = mg_node_alloc(mi);
    d=(double*)malloc(sizeof(double));
    *d=func.evaluate(mi, func.data);
    node->info=d;
    mpa_monotone_graph_add_node(omg,node);
    
    /*####### iterate and add nodes one by one to the monotone graph #######*/
    mpa_range focus=mpa_range_make(1,0);
    for (int j=0; j<capacity-1; j++) {
        //printf("j======%d\n",j);
        node = omg->nodes[j];
        node->IN = 1;
        mi = node->miPtr;
   
        for (int k=0; k<dim; k++) {
            int  is_admissible = mg_node_successor_all_predecessors_in(node, k);
            if (is_admissible) {
                succ_mi   = mpa_multi_index_plus_ei(mi,k);
                succ_node = mg_node_alloc(succ_mi);
                d=(double*)malloc(sizeof(double));
                *d= func.evaluate(succ_mi,func.data);
                //node->info = (double*) malloc(sizeof(double));
                succ_node->info=d;
                mg_node_set_all_predecessors(succ_node, node, k);
                mpa_monotone_graph_insert_node(omg, focus, succ_node,mg_node_compare_suite);
                focus.length++;
            }
        }
        focus.position++;
        focus.length=omg->size-focus.position;
    }
    //	mpa_monotone_graph_resize(omg, capacity);
    
    /*
    FILE * enf = fopen("/Users/abdellahchkifa/Documents/Programmation/programmationC/
     bulk_chase/cas 16 psi/results/rhoFileME.txt", "w");
    
    	for (int i=0; i<omg->size; i++){
    		double rhos=*((double*)omg->nodes[i]->info);
    		fprintf(enf,"%e\n",0.00158084/0.05*rhos);
    		free((double*)omg->nodes[i]->info);
    	}*/
    return omg;
}
