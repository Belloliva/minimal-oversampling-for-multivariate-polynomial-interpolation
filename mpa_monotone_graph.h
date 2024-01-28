/*
 *  mpa_monotone_graph.h
 *  mpa
 *
 *  Created by Abdellah CHKIFA on 24/02/11.
 *  Copyright 2011 oliva Electronics. All rights reserved.
 *
 */
#ifndef __mpa_monotone_graph_h__
#define __mpa_monotone_graph_h__

#include "mpa_multi_index.h"
#include "mpa_function.h"
#include "mpa_stdlib.h"
#include <stdarg.h>

/*! \struct mg_node
 \brief the monotone graph node description.
 
 describe the node in a monotone graph, each node has the member:
 info, a void pointer to hold the information used in an algorithm
 miPtr, multi-index pointer
 successors, array of node which are successors (or childs) of the current node
 predecessors, array of node which are successors (or parents) of the current node
 */

typedef struct mg_node mg_node;

struct mg_node{
    mpa_multi_index* miPtr;
    size_t mi_abs;
    mg_node** successors;
    mg_node** predecessors;
    int birthIteration;
    int IN;
    void *info;
};

/*! \struct mpa_monotone_graph
 \brief the monotone graph node description.
 
 a dynamic array of the graph nodes:
 capacity, the capacity the the dynamic array.
 size, the size of the dynamic array, or the count of nodes that the array hold
 nodes, the array holding the nodes
 */

typedef struct{
    size_t capacity;
    size_t size;
    mg_node** nodes;
}mpa_monotone_graph;



#pragma mark -
#pragma mark Monotone graph nodes functions

/*!
 * Creates a node, returning a pointer to a newly initialized node struct. An array of nodes pointer of length the multi-index length is allocated, and stored in the successors component and an array of nodes pointer of length the multi-index size is allocated, and stored in the predecessors component . These arrays are “owned” by the node, and will be deallocated when the node is deallocated. info is equal to NULL.
 
 * @param miPtr the multi-index pointer
 * @return pointer the new allocated node
 */
mg_node* mg_node_alloc(mpa_multi_index* miPtr);

/*!
 *Free a previously allocated mg_node node.
 * @param node the mg_node to be freed.
 */
void mg_node_free(mg_node* node);

/*!
 Compare the two nodes pointed to by comparing corresponding multi-indices absolute values. The result is >0 if node1 multi-index
 abslolute value is greater than node2 multi-index abslolute value, 0 is equality and < 0 else.
 The node pointer are casted to void because this function is used with sorting algorithm.
 
 *@param node1 the first node.
 *@param node2 the second node.
 *@return and int.
 */
int mg_node_compare_abs            (const void* node1,const void* node2);



/*!
 Set succ as the \f$i^{th}\f$ successors or "child" of node,
 *@param node the first node.
 *@param succ the successor node.
 *@param i the successor order.
 */
void mg_node_set_successor(mg_node* node, mg_node* succ, int i);

/*!
 Set pred as the \f$p^{th}\f$ predeccessors or "parent" of node,
 *@param node the first node.
 *@param pred the predecessor node.
 *@param p the predecessor order.
 */
void mg_node_set_predecessor(mg_node* node, mg_node* pred, int p);

/*!
 Set all the predecessors of a node by only knowing one predecessor, this is always valid in a monotone graph.
 *@param node the first node.
 *@param aPred a predecessor node.
 */
void mg_node_set_all_predecessors(mg_node* node, mg_node* aPred, size_t p);


/*!
 Check that a node has all its predecessor already in the monotone graph.
 
 *@param node the node.
 *@return the boolean value.
 */
int  mg_node_all_predecessors_in(const mg_node* node);
int  mg_node_successor_all_predecessors_in(const mg_node* node, size_t k);


//void mg_node_printf(mg_node* node);
//mg_node* mg_node_find_successor(mg_node* node,size_t k,mg_node** margin,size_t size);
//mg_node* mg_node_find_predecessor(mg_node* node,size_t k,mg_node** margin,size_t size);
//mpa_search_pair  mg_node_insert(mg_node* node, mg_node** nodeList,size_t size);

#pragma mark -
#pragma mark Monotone graph functions

/*!
 * Creates a monotone graph, returning a pointer to a newly initialized monotone graph struct. An array of nodes pointer of length the
 * capacity length is allocated, and stored in are the nodes component. This array is “owned” by the graph, and will be deallocated
 * when the node is deallocated.
 
 * @param capacity the initial capacity of the nodes dynamic array
 * @return pointer the new allocated monotone graph
 */
mpa_monotone_graph* mpa_monotone_graph_alloc(size_t capacity);

/*!
 * Free a previously allocated monotone graph. All the nodes of the graph are freed along the way.
 * @param omg the monotone graph to be freed.
 */
void mpa_monotone_graph_free(mpa_monotone_graph* omg);


/*!
 * Creates a monotone graph , returning a pointer to a newly initialized monotone graph struct. The graph
 is constructed as the nodes that correspond to the largest value of the function evaluated at the node
 * @param dim the common length of all the multi-indices of the node
 * @param capacity the capacity and the size of the constructed graph
 * @param func a function to be evaluated at the node
 * @return pointer the new allocated monotone graph
 */
mpa_monotone_graph* mpa_monotone_graph_envellope(size_t dim,
                                                 size_t capacity,
                                                 mpa_function func);


/*!
 * Resize a graph,if the new capacity is less than the old capacity, the nodes beyond new capacity are freed, if new capacity is
 * greater than the old capacity, the array of nodes is reallocared and all the nodes copied to it.
 * @param omg the monotone graph
 * @param newCapacity the new capacity graph
 */
void mpa_monotone_graph_resize(mpa_monotone_graph* omg,size_t newCapacity);

/*!
 * Add a node to the graph, if the graph is full, it is resized to double the current capacity. The graph size is incremend by one
 * @param omg the monotone graph
 * @param node the  node
 */
void mpa_monotone_graph_add_node(mpa_monotone_graph* omg,mg_node* node);

/*!
 * Add a node to the graph in a range, we use the function cmp to check where to include the node. If the graph is full, it is resized
 * to double the current capacity. The graph size is incremend by one
 * @param omg the monotone graph
 * @param node the  node to be inserted
 * @param r the range where to insert the node
 * @param cmp the comparaison predicate
 */
void mpa_monotone_graph_insert_node(mpa_monotone_graph* omg,mpa_range r, mg_node* node,
                                    int (*cmp)(const void *, const void *));




#endif /*__mpa_monotone_graph_h__*/
