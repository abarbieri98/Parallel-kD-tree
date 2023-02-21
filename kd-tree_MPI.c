#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "mpi.h"
#include "util_functions.h"

#define NDIM 2

/*
	TODO:
		- clean code and functions from extra variables
		- move shared functions with kd-tree.c to sorting.h and sorting.c
*/

int build_tree(kpoint* dataset, struct kdnode_mpi* tree, int* indexes, int n_points, int node_index, int depth, int size, int rank);
void extract_points(kpoint* dataset, kpoint* arr, int* indexes, int n_points);
int get_sender(int me);
void init_points(kpoint* arr, int n);
void init_nodes(struct kdnode_mpi* arr, int n);
void join_splits(struct kdnode_mpi* big_arr, struct kdnode_mpi* small_arr, int size_big, int size_small);


int main(int argc, char **argv){

    if(argc != 2){
        printf("Error: invalid number of command line args.\nUsage: [n_points]\n");
        return 1;
    }

    int n_points = atoi(argv[1]);
    
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Status status;

    MPI_Datatype mpi_kpoint;
    MPI_Type_contiguous(2,MPI_DOUBLE,&mpi_kpoint);
    MPI_Type_commit(&mpi_kpoint);

    kpoint *dataset;
    kpoint *local_data;

    // Here we define a custom MPI type to send 
    // nodes

	MPI_Datatype mpi_kdnode_mpi;
	int lengths_node[4] = {1,1,1,1};
	MPI_Aint displacements_node[4];
	
	struct kdnode_mpi node;

	MPI_Aint base_address;
	MPI_Get_address(&node, &base_address);
	MPI_Get_address(&node.split.coord[0], &displacements_node[0]);
	MPI_Get_address(&node.axis, &displacements_node[1]);
	MPI_Get_address(&node.left, &displacements_node[2]);
	MPI_Get_address(&node.right, &displacements_node[3]);

	displacements_node[0] = MPI_Aint_diff(displacements_node[0], base_address);
	displacements_node[1] = MPI_Aint_diff(displacements_node[1], base_address);
	displacements_node[2] = MPI_Aint_diff(displacements_node[2], base_address);
	displacements_node[3] = MPI_Aint_diff(displacements_node[3], base_address);
	
	MPI_Datatype types_node[4] = {mpi_kpoint, MPI_INT, MPI_INT, MPI_INT };
	MPI_Type_create_struct(4, lengths_node, displacements_node, types_node, &mpi_kdnode_mpi);
	MPI_Type_commit(&mpi_kdnode_mpi);
    
    if(rank==0){
        
        //  1st step: data allocation and reading
        dataset = (kpoint*) malloc(sizeof(kpoint)*n_points);
        init_points(dataset,n_points);
        read_dataset(dataset, n_points);

        //  2nd step: allocate variables to be used
        //  in the next steps
        struct kdnode_mpi* local_tree = (struct kdnode_mpi*) malloc(sizeof(struct kdnode_mpi)*n_points);
		init_nodes(local_tree, n_points);
		int* indexes = (int*) malloc(sizeof(int)*n_points);
        create_indexes(n_points,indexes);
        double start = MPI_Wtime();

        //  3rd step: Tree building and data transfer
        //  to other processes
        build_tree(dataset,local_tree,indexes,n_points,0,0,size,rank);

        //  4th step: gather nodes from the other processes and join them
        //  to the local tree
        int n_iter =log2(size);
        
            for(int iter = 0; iter < n_iter; iter++){
                int n_arriving = 0;
                int sender = rank + pow(2,n_iter-iter-1);
                MPI_Recv(&n_arriving, 1, MPI_INT, sender, sender, MPI_COMM_WORLD, &status);
                
                struct kdnode_mpi* temp_storage = (struct kdnode_mpi*) malloc(sizeof(struct kdnode_mpi)*n_arriving);
                MPI_Recv(temp_storage,n_arriving, mpi_kdnode_mpi, sender,sender, MPI_COMM_WORLD, &status);
                join_splits(local_tree,temp_storage,n_points,n_arriving);
                free(temp_storage);
            }
        double end = MPI_Wtime();
        printf("Elapsed time to build tree: %10.5lf\n", end - start);
    }
    else{

        // 1st step: receive the total number of points and the points
        //  themselves 
        int* informations = (int*) malloc(sizeof(int)*2);
        int sender = get_sender(rank);
        MPI_Recv(informations,2,MPI_INT,sender, sender, MPI_COMM_WORLD,&status);
        kpoint* local_data = (kpoint*) malloc(sizeof(kpoint)*informations[0]);
        init_points(local_data,informations[0]);
        MPI_Recv(local_data,informations[0],mpi_kpoint,sender,sender, MPI_COMM_WORLD, &status);

        //  2nd step: variables allocation for the later steps
        struct kdnode_mpi* local_tree = (struct kdnode_mpi*) malloc(sizeof(struct kdnode_mpi)*informations[0]);
		init_nodes(local_tree, informations[0]);
        int* indexes = (int*) malloc(sizeof(int)*informations[0]);
        create_indexes(informations[0],indexes);

        // 3rd step: Build the local tree 
        build_tree(local_data,local_tree,indexes,informations[0],0,informations[1],size,rank);

        // 4th step: depending on the process rank, send or receive nodes
        //  from other processes
        if(rank >= size/2){
            int receiver = get_sender(rank);
            MPI_Send(&informations[0],1,MPI_INT, receiver, rank, MPI_COMM_WORLD);
            MPI_Send(local_tree, informations[0], mpi_kdnode_mpi, receiver , rank, MPI_COMM_WORLD);
        }
        else{

            int n_iter = log2(size)-floor(log2(rank))-1 ;

            for(int iter = 0; iter < n_iter; iter++){
                int n_arriving = 0;
                int sender = rank + pow(2,log2(size)-iter-1);
                MPI_Recv(&n_arriving, 1, MPI_INT, sender, sender, MPI_COMM_WORLD, &status);
                struct kdnode_mpi* temp_storage = (struct kdnode_mpi*) malloc(sizeof(struct kdnode_mpi)*n_arriving);
                MPI_Recv(temp_storage,n_arriving, mpi_kdnode_mpi, sender,sender, MPI_COMM_WORLD, &status);
                join_splits(local_tree,temp_storage,informations[0],n_arriving);
                free(temp_storage);
            }

            int receiver = get_sender(rank);
            MPI_Send(&informations[0],1,MPI_INT, receiver, rank, MPI_COMM_WORLD);
            MPI_Send(local_tree, informations[0], mpi_kdnode_mpi, receiver , rank, MPI_COMM_WORLD);
        }
        
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}

int build_tree(kpoint* dataset, struct kdnode_mpi* tree, int* indexes, int n_points, int node_index, int depth, int size, int rank){

    /*
        Function for tree-building inside each processor.
        The function needs an array of nodes, where it will
        store the computed nodes in a post-order fashion.

        Summarizing, the steps are:
            1 - choose axis
            2 - find splitting point
            3 - split left and right points
            4 - accordingly to the depth, send data points
                or go on with the creation of nodes
    */

    //  1st step: choose axis
    int my_axis = (depth)%NDIM; 
    struct kdnode_mpi* this_node = (struct kdnode_mpi*) malloc(sizeof(struct kdnode_mpi));
    
    int last_node;
    tree[node_index].axis = my_axis;

    //  Terminal node condition
    if(n_points == 1){
        
		tree[node_index].split = dataset[indexes[0]];
		tree[node_index].left = -1;
		tree[node_index].right = -1;
        free(indexes);
        last_node = node_index;
    }
    else{
        
        // 2nd step: find splitting point 
        tree[node_index].split = best_split(dataset, indexes, n_points/2, my_axis, 0, n_points -1);

        // 3rd step: split left and right points for the next nodes
        n_points--;
        int nr = (n_points)/2;
        int nl = (n_points%2 == 0? nr : nr+1);
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
        
        //  4th step: accordingly to the depth, send data points
        //  or go on with the creation of nodes

        if(pow(2,depth) < size){

            //  If the condition is true, right points will be sent
            //  to another process while the sender process will
            //  go on with the left ones

            MPI_Datatype mpi_kpoint;
            MPI_Type_contiguous(2,MPI_DOUBLE,&mpi_kpoint);
            MPI_Type_commit(&mpi_kpoint);
            

            kpoint* right_points = (kpoint*) malloc(sizeof(kpoint)*nr);
            init_points(right_points,nr);
            extract_points(dataset, right_points, right_indexes, nr);
            int receiver = rank + pow(2,depth);
            int* info = (int*) malloc(sizeof(int)*2);
            
            info[0] = nr;
            info[1] = depth+1;
            MPI_Send(info,2,MPI_INT, receiver, rank, MPI_COMM_WORLD);
            free(info);
            MPI_Send(right_points, nr, mpi_kpoint, receiver, rank, MPI_COMM_WORLD);
            
			tree[node_index].left = node_index + 1;
            build_tree(dataset, tree, left_indexes, nl, node_index+1, depth+1, size, rank);

        }
        else{

            //  If the process does not have to send other points
            //  it will go on building the nodes of the tree, iteratively calling
            //  the build_tree function

			tree[node_index].left = node_index + 1;
            last_node = build_tree(dataset,tree, left_indexes, nl, node_index+1, depth+1, size, rank);
            
            if(nr>0){
				tree[node_index].right = last_node + 1;
                last_node = build_tree(dataset,tree, right_indexes, nr, last_node+1, depth+1, size, rank);
            }
			else{
				tree[node_index].right = -1;
			}
        }
    }   
        
    return last_node;
}

void extract_points(kpoint* dataset, kpoint* arr, int* indexes, int n_points){
    /*
        Create needed dataset "arr" from a list of indexes
    */
    for(int i = 0; i<n_points; i++)
        arr[i] = dataset[indexes[i]];

}

int select_axis(int size){
    int res = ((int)log2(size))%2;
    res = (res == 0 ? 1 : 0);
    
    return res;
}

int get_sender(int me){

    //  Given a process return the process that
    //  will send the computed nodes in the 4th 
    //  step of the main functon
    return me - pow(2,floor(log2(me)));
}

void init_points(kpoint* arr, int n){

    //Initialize an array of kpoints

    for(int i = 0; i < n; i++){
        arr[i].coord[0] = -1;
		arr[i].coord[1] = -1;
    }
}

void init_nodes (struct kdnode_mpi* arr, int n){

    //initialize an array on kdnodes

	for(int i = 0; i < n; i++){
        arr[i].axis = -1;
    }
}

void join_splits(struct kdnode_mpi* big_arr, struct kdnode_mpi* small_arr, int size_big, int size_small){

    //  Function to join two array of kdnondes, used in the 
    //  4th step of the main function
    
    int diff = 0;
    while(big_arr[diff].axis != -1){
        diff++;
    }
	big_arr[0].right = diff;

    for(int i =0; i < size_small; i++){
        big_arr[i+diff] = small_arr[i];
		if(big_arr[i+diff].left > 0){
			big_arr[i+diff].left += size_big; 
		}
		if(big_arr[i+diff].right > 0){
			big_arr[i+diff].right += size_big; 
		}
    }
}