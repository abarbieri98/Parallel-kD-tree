#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<omp.h>
#include<time.h>
#include "util_functions.h"

#define NDIM 2

struct kdnode* build_tree (kpoint* dataset, int* indexes, int n_points, int depth );

int main(int argc, char **argv){

    if(argc != 3){
        printf("Error: invalid number of command line args.\nUsage: [n_points] [n_threads]\n");
        return 1;
    }

    int n_points = atoi(argv[1]);
    int p = atoi(argv[2]); //Number of threads to use
    #if defined(_OPENMP)
    omp_set_dynamic(0);
    omp_set_num_threads(p);
    #endif

    kpoint* dataset = (kpoint*) malloc(sizeof(kpoint)*n_points);
    read_dataset(dataset, n_points);

    int* indexes = malloc(sizeof(int)*n_points);
    create_indexes(n_points,indexes);

    struct kdnode *tree = (struct kdnode*) malloc(sizeof(struct kdnode));
    
    #if defined(_OPENMP)
    double start, end; 
    start = omp_get_wtime();
    #pragma omp parallel shared(dataset)
        {
            #pragma omp single 
            {
                tree = build_tree(dataset, indexes, n_points, 0);
            }
        }
    end = omp_get_wtime();
    printf("Elapsed time to build tree: %10.5lf\n",(end-start));
    #else
    clock_t start,end;
    start = clock();
    tree =  build_tree(dataset, indexes, n_points, 0);
    end = clock();
    printf("Elapsed time to build tree: %10.5lf\n",(double) (end-start)/CLOCKS_PER_SEC);
    #endif
    
    free(tree);
    free(dataset);
    return 0;
}

struct kdnode* build_tree (kpoint* dataset, int* indexes, int n_points, int depth){

    /*
        Builds a node of the tree.
        Sequentially, it will:
        - allocate the space for the node 
        - choose the splitting axis 
        - search the best splitting point
        - assign the observations for its children building
        - task the creation of its children

        For 1-point nodes the last two steps
        will be ignored.
        
    */
    
    struct kdnode* this_node = (struct kdnode*) malloc(sizeof(struct kdnode));
    int my_axis = (depth+1)%NDIM; 
    this_node->axis = my_axis;
    
    
    if(n_points == 1){
        kpoint split;
        split = dataset[indexes[0]];
        this_node->split = split;
        this_node->left = NULL;
        this_node->right = NULL;
        free(indexes);
    }
    else{

        //sort_index(dataset, indexes, 0, n_points-1, axis);
        
        kpoint split; 
        int med_ind = n_points/2;
        split = best_split(dataset, indexes, med_ind, my_axis, 0, n_points-1);
        this_node->split = split;
        n_points--;
        int nr = (n_points)/2;
        int nl = nr;
        if (n_points%2 != 0){
            nl++;
        }
        int* left_indexes;
        int* right_indexes;

        if(nr > 0){
            left_indexes = (int*) malloc(sizeof(int)*nl);
            right_indexes = (int*) malloc(sizeof(int)*nr);

            assign_lr(indexes, left_indexes, right_indexes, nl, nr);
        }
        else{
            left_indexes = (int*) malloc(sizeof(int)*nl);
            left_indexes[0] = indexes[0];
        }
        free(indexes);
            //int me = omp_get_thread_num();
        #if defined(_OPENMP)
        {
                if(nr > 0){
                    #pragma omp task shared(this_node, dataset) 
                    {
                        this_node->right = build_tree(dataset,right_indexes,nr, depth+1);
                    }
                }
                else{
                    this_node->right = NULL;
                }

                this_node->left = build_tree(dataset,left_indexes,nl, depth+1);
        }
        
        #else
        
        this_node->left = build_tree(dataset,left_indexes,nl, depth+1);
        if(nr >0){
            this_node->right = build_tree(dataset,right_indexes,nr, depth+1);
            
        }
        else{
            this_node->right = NULL;
        }
        #endif
    }
    
    return this_node;
}  
