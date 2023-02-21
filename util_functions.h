#define NDIM 2
typedef struct {double coord[NDIM];}kpoint;

struct kdnode{
int axis;// the splitting dimension
kpoint split;// the splitting element
struct kdnode *left, *right;// the left and right sub-trees
};

struct kdnode_mpi{
    int axis;// the splitting dimension
    kpoint split;// the splitting element
    int left,right;// the left and right sub-trees
};
void sort_index (kpoint *array, int* indexes, double left, double right, int axis);
void merge (kpoint* arr, int* indexes, int l, int m, int r, int axis);
void print_points(kpoint* arr, int n);
void create_indexes(int max, int* indexes);
int partition(kpoint *arr, int* indexes, int axis, int l, int r, int p_index);
void sort_points(kpoint *arr, int low, int high, int axis);
void print_points_x(kpoint* arr, int n, int* indexes);
void read_dataset (kpoint *dataset, int n_points);
kpoint best_split (kpoint* dataset, int* indexes, int median_index, int axis, int l, int r);
void assign_lr(int* indexes, int* l, int* r, int nl, int nr);