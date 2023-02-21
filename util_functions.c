#include<stdlib.h>
#include<stdio.h>
#include<omp.h>
#include "util_functions.h"

void swap(int a, int b, int *arr)
{
    int t = arr[a];
    arr[a] = arr[b];
    arr[b] = t;
}

int partition(kpoint *arr, int* indexes, int axis, int l, int r, int p_index){
    
    double pivot = arr[indexes[p_index]].coord[axis];
    swap(p_index,r, indexes);
    p_index = l;

    for(int i=l; i<r; i++){
        if(arr[indexes[i]].coord[axis] <= pivot){
            swap(i, p_index, indexes);
            p_index++;
        }
    }

    swap(p_index,r,indexes);
    return p_index;
}

kpoint best_split (kpoint* dataset, int* indexes, int median_index, int axis, int l, int r){

    // Choose the median point of an array
    // (Note: in case of even number of observations
    //  it will return the (n+2)/2 -th point in the dataset)

    if(l==r){
        return(dataset[indexes[l]]);
    } 
    int p_index = (l + r) / 2;

    p_index = partition(dataset, indexes, axis, l, r, p_index);
    
    if(median_index == p_index){
        return dataset[indexes[p_index]];
    }
    else if (median_index < p_index) {
        return best_split(dataset, indexes, median_index, axis ,l, p_index-1);
    }
    else{
        return best_split(dataset, indexes, median_index, axis, p_index+1, r);
    }

}



void print_points(kpoint* arr, int n){
    for(int i = 0; i < n ; i++){
        printf("(%lf,%lf)\n", arr[i].coord[0], arr[i].coord[1]);
    }
    printf("\n");
}

void print_points_x (kpoint* arr, int n, int* indexes){
    for(int i = 0; i < n ; i++){
        printf("(%lf,%lf)\n", arr[indexes[i]].coord[0], arr[indexes[i]].coord[1]);
    }
    printf("\n");
}

void create_indexes(int max, int* indexes){
    for(int i = 0; i < max; i++){
        indexes[i]=i;
    }
}

void read_dataset (kpoint *dataset, int n_points){
    
    // simple function to read a csv file and save points
    // to an array
    FILE* f = fopen("./dataset.csv", "r");

    for(int i = 0; i < n_points; i++){
        fscanf(f,"%lf,%lf", &dataset[i].coord[0], &dataset[i].coord[1]);
        
    }
    
}


void assign_lr(int* indexes, int* l, int* r, int nl, int nr){

    //  Function to assign the left and right indexes
    //  of the points, defining the points of the two
    //  children nodes of a node

    for(int i = 0; i < nr; i++){
        l[i] = indexes[i];
        r[i] = indexes[nl+i+1];
    }
    if (nl > nr){
        l[nl-1] = indexes[nl-1];
    }
}
