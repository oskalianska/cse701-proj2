#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <stdlib.h>

// All possible values node types:
//
// Node that was never intialized and represents an invalid node that must be
// ignored if encountered in any of the calculations
#define NODE_NOT_IN_USE 0
// Node type that holds data; these types of nodes are combined into clusters
#define NODE_LEAF 1
// Node that is a non-data element. These nodes are used to combine two leaf nodes
// or other merged nodes into one
#define NODE_MERGED_ROOTS 2

// Linkage critereas (the way of calculation of distances between clusters):
//
// Use computed average value between the clusters
#define LINKAGE_AVG 'a'
// Use computed distance between the centeroids of the clusters
#define LINKAGE_CENTROID 't'
// Use maximum distance between the classes as a linkage criteria
#define LINKAGE_COMPLETE 'c'
// Use minimum distance between the classes as a linkage criteria
#define LINKAGE_SINGLE 's'

// Pre-processor macro, a shortcut for calloc without the need to specify the size of
// each element in the sequence
#define allocate_memory(SIZE, T) (T *)calloc(SIZE, sizeof(T))

// If PRINT_HIERARCHY is defined, then hierarchy information will be sent to standard
// output when algorithm finishes execution
#define PRINT_HIERARCHY
// If PRINT_CLUSTERS is defined, then all information about identified clusters and all
// elements, contained withing each of the clusters, will be sent to standard output
// when algorithm finishes execution
#define PRINT_CLUSTERS

// Constant to define maximum length of a label of a leaf element (data-element)
const int MAX_LABEL_SIZE = 16;

// Define shortcut to refer to struct types without specifying "strcut" keyword
typedef struct cluster_struct cluster_type;
typedef struct node_struct node_type;
typedef struct neighbour_struct neighbour_type;
typedef struct labeled_data_struct labeled_data_type;

// Define a signature for calculation of the distance between two nodes
// Based on required method of distance calculation in command line arguments, a specific
// actual implementaiton of the function will be used
float (*dist_float_ptr)(float **, const int *, const int *, int, int);

// Definition of a cluster type
struct cluster_struct
{
        int items_count;
        int clusters_count;
        int nodes_count;
        node_type *nodes;
        float **dist;
};

// Structure to hold information about data values of the leaf node
typedef struct data_struct
{
        float x, y;
} coord_type;

// Data structure to hold full information about a node
struct node_struct
{
        int type;                   // type of node
        int is_root;                // marker whether node is a root element in the cluster
        int height;                 // height of the node in the cluster
        coord_type centroid;        // calculated values for a center of the node
        char *tag;                  // label of the node; each data value has its own unique label
        int *merged;                // holds indexes to the two merged values
        int items_count;            // count of element inside of the node
        int *items;                 // addresses of all leaf elements, belonging to the node
        neighbour_type *neighbours; // pointers to the neighbour nodes of the current node
};

// Data type to hold information about neighboours of a node
struct neighbour_struct
{
        int point;                       // index of the neighbouring node
        float dist;                      // distance to neighbour's node
        neighbour_type *previous, *next; // pointers to the neighbour nodes
};

// Type that describe poinr's coordinates and it's label
struct labeled_data_struct
{
        coord_type coordinate;    // loading a coordinate data
        char tag[MAX_LABEL_SIZE]; // adding label of data
};

// Function to calculate the Euclidean distance between two coordinates
float calculate_dist(const coord_type *a, const coord_type *b)
{
        return sqrt(pow(a->x - b->x, 2) + pow(a->y - b->y, 2));
}

// Function to calculate and populate the matrix of distances 'm' between the data elements in 'items' array
void fill_distances(float **m, int items_count, const labeled_data_type items[])
{
        for (int i = 0; i < items_count; ++i)
        {
                for (int j = 0; j < items_count; ++j)
                {
                        m[i][j] = calculate_dist(&(items[i].coordinate), &(items[j].coordinate));
                        m[j][i] = m[i][j];
                }
        }
}

// Create a matrix of distances for all data elements in the array 'items', based on the method
// for calculating distances, specified as a command line argument
float **create_dist_matrix(int items_count, const labeled_data_type items[])
{
        // allocate memory for the pointers to each row of the matrix
        float **m = allocate_memory(items_count, float *);
        if (m)
        {
                for (int i = 0; i < items_count; ++i)
                {
                        // allocate memoty for the row number 'i' in the resulting array
                        m[i] = allocate_memory(items_count, float);
                        if (!m[i])
                        {
                                // if we run our of memory or exceed the limit, available to the appliation,
                                // show an error the the use and clean-up already allocated memory before
                                // existing the for loop
                                fprintf(stderr, "Failed to allocate memory for distance matrix row.\n");
                                items_count = i;
                                for (i = 0; i < items_count; ++i)
                                {
                                        free(m[i]);
                                }
                                free(m);
                                m = NULL;
                                break;
                        }
                }
                // memory allocation was successfull
                if (m)
                {
                        // fill our matrix with distances between elements in 'items' array
                        fill_distances(m, items_count, items);
                }
        }
        else
        {
                // print error to stanrad error output
                fprintf(stderr, "Failed to allocate memory for distance matrix.\n");
        }
        return m;
}

// Function to find minimal distance between any two elements if 'a' and 'b'
float single_link(float **dist, const int a[], const int b[], int m, int n)
{
        // start with largest possible float value
        // distance between any two nodes must be smaller than this value
        // and must override it on the first pass
        float res = FLT_MAX, d;
        for (int i = 0; i < m; ++i)
        {
                for (int j = 0; j < n; ++j)
                {
                        d = dist[a[i]][b[j]];
                        if (d < res)
                                res = d;
                }
        }
        return res;
}

// Function to find the maximum distance between any two elements if 'a' and 'b'
float complete_link(float **dist, const int a[], const int b[], int m, int n)
{
        // start with MAX distance 0
        // distance between any two nodes must be larger or equal to this value
        // and ,if larger, it must override it
        float d, res = 0.0;
        for (int i = 0; i < m; ++i)
                for (int j = 0; j < n; ++j)
                {
                        d = dist[a[i]][b[j]];
                        if (d > res)
                                res = d;
                }
        return res;
}

// Function to find the average of distance between 'a' and 'b' based on all the elements in 'a' and 'b'
float avg_link(float **dist, const int a[], const int b[], int m, int n)
{
        float res = 0.0;
        for (int i = 0; i < m; ++i)
                for (int j = 0; j < n; ++j)
                        res += dist[a[i]][b[j]];
        return res / (m * n);
}

// Dummy function to find the centroid distance between 'a' and 'b'. It is en empty function 
// as centroid will be calculated by its own 'calculate_dist' funciton
float centr_link(float **dist, const int a[], const int b[], int m, int n)
{
        return 0;
}

// Get distance between two nodes based on the specified linkage criteria
float get_dist(cluster_type *c, int index, int point)
{
        // for two leaf elements get distance from the cluster's matrix of distances
        if (index < c->items_count && point < c->items_count) {
                return c->dist[index][point];
        }
        // for the elements, where at least one is not-leaf one, distance must be calculated
        else
        {
                node_type *a = &(c->nodes[index]);
                node_type *b = &(c->nodes[point]);
                // centroid linkage between two nodes required a special calculation
                if (dist_float_ptr == centr_link)
                {
                        return calculate_dist(&(a->centroid), &(b->centroid));
                }
                // for any other type of linkage, use its appropriate function
                else
                {
                        return dist_float_ptr(c->dist, a->items, b->items, a->items_count, b->items_count);
                }
        }
}

// Function to free up memory, used by a neighbour element.
// It will free up memory for the element as well as all linked elements.
void free_memory(neighbour_type *n)
{
        neighbour_type *t;
        while (n)
        {
                t = n->next;
                free(n);
                n = t;
        }
}

// Function to free up memory, taken up by all nodes of the cluster.
void free_nodes_memory(cluster_type *c)
{
        for (int i = 0; i < c->nodes_count; ++i)
        {
                node_type *node = &(c->nodes[i]);
                // if node has text label, free up its memory for the label
                if (node->tag)
                        free(node->tag);
                // if node has merged nodes, free up all memory taken by them
                if (node->merged)
                        free(node->merged);
                // if node has any leaf elements, belonging to it, free up their memory
                if (node->items)
                        free(node->items);
                // if node has any neighbouring nodes, free up their memory
                if (node->neighbours)
                        free_memory(node->neighbours);
        }
        free(c->nodes);
}

// Function to free up memory, taken by the cluster. It will free up distance matrix 
// as well as all nodes in the cluster
void free_cluster_memory(cluster_type *c)
{
        // if c is not a NULL-pointer, proceed with the clean up
        if (c)
        {
                // if nodes is not a NULL pointer, free up nodes memory
                if (c->nodes)
                        free_nodes_memory(c);
                // if distance matrix is not a NULL pointer, clean it up
                if (c->dist)
                {
                        // clean up memory taken up by each row
                        for (int i = 0; i < c->items_count; ++i)
                                free(c->dist[i]);
                        // free up the list of pointers to each row
                        free(c->dist);
                }
                // free up memory, taken up by the cluster
                free(c);
        }
}

// Associate node with neigbours, as a node following them
void ins_before(neighbour_type *curr, neighbour_type *neighbours, node_type *n)
{
        // set next element to the current one
        // is neighbours are inserted before curr, it means that next for neighbours 
        // will be current element 
        neighbours->next = curr;

        if (curr->previous)
        {
                // if current element has non-empty previous element, replace it with neighbours
                // and fix relation so that old previous element becomes second previous now
                curr->previous->next = neighbours;
                neighbours->previous = curr->previous;
        }
        else
        {
                // if current element has an empty previous element, then set node neighbours
                n->neighbours = neighbours;
        }

        // associcate 'curr' element with the 'neighbours' as its previous element
        curr->previous = neighbours;
}

// Associate node with neigbours, as a node preceding them
void ins_after(neighbour_type *curr, neighbour_type *neighbours)
{
        neighbours->previous = curr;
        curr->next = neighbours;
}

// Associate neighbours, based on distance between them
void ins_sorted(node_type *node, neighbour_type *neighbours)
{
        neighbour_type *temp = node->neighbours;
        // by checking next neighbours of the node, if we find node with larger distance than passed neighbours,
        // then insert passed neighbours before that node
        while (temp->next)
        {
                // check if distance is larger or same, than passed neighbours
                if (temp->dist >= neighbours->dist)
                {
                        ins_before(temp, neighbours, node);
                        // stop execution when smaller or equal distance found
                        return;
                }
                // set element for the next iteration
                temp = temp->next;
        }
        // compare distance to last element
        if (neighbours->dist < temp->dist)
        {
                // insert neighbour as precedding last element
                ins_before(temp, neighbours, node);
        }
        else
        {
                // insert neighbour as following last element
                ins_after(temp, neighbours);
        }
}

// Add a node with the specified index to a cluster, placing it in the sorted order 
neighbour_type *add_neighbour(cluster_type *cluster, int index, int point)
{
        neighbour_type *neighbour = allocate_memory(1, neighbour_type);
        if (neighbour)
        {
                neighbour->point = point;
                neighbour->dist = get_dist(cluster, index, point);
                node_type *node = &(cluster->nodes[index]);
                if (node->neighbours)
                        // if node has neoughbours, insert in the sorted order
                        ins_sorted(node, neighbour);
                else
                        // if it's the first node in cluster, add neighbours as is
                        node->neighbours = neighbour;
        }
        else
        {
                // print error to stanrad error output
                fprintf(stderr, "Failed to allocate memory for neighbour node.\n");
        }

        return neighbour;
}

// Function processes node with index 'i' in the cluster
cluster_type *update_neighbour_nodes(cluster_type *c, int i)
{
        node_type *node = &(c->nodes[i]);
        // node 'i' must be a part of the clsuter
        if (node->type == NODE_NOT_IN_USE)
        {
                // if it is not a part of the cluster, print an error message
                fprintf(stderr, "Invalid cluster node at index %d.\n", i);
                c = NULL;
        }
        else
        {
                int root_clusters_seen = 1, point = i;
                while (root_clusters_seen < c->clusters_count)
                {
                        node_type *tmp = &(c->nodes[--point]);
                        if (tmp->type == NODE_NOT_IN_USE)
                        {
                                // print an error message if node is not a part of the cluster
                                fprintf(stderr, "Invalid cluster node at index %d.\n", i);
                                c = NULL;
                                break;
                        }
                        if (tmp->is_root)
                        {
                                ++root_clusters_seen;
                                add_neighbour(c, i, point);
                        }
                }
        }
        return c;
}

// Pre-processor macro to initialize a leaf node 'n' from cluster 'n'
#define initialize_leaf(c, n, item, label_size)         \
        {                                               \
                strncpy(n->tag, item->tag, label_size); \
                n->centroid = item->coordinate;         \
                n->type = NODE_LEAF;                    \
                n->is_root = 1;                         \
                n->height = 0;                          \
                n->items_count = 1;                     \
                n->items[0] = c->nodes_count++;         \
        } 

// Function to add leaf 'item' to the cluster 'c'.
// It returns newly added leaf node.
node_type *add_leaf(cluster_type *c, const labeled_data_type *item)
{
        // get pointer to the last element in the cluster
        node_type *leaf = &(c->nodes[c->nodes_count]);
        // allocate memory for the tag
        int length = strlen(item->tag) + 1;
        leaf->tag = allocate_memory(length, char);
        // if memory was successfully allocated
        if (leaf->tag)
        {
                // proceed with allocating memory for the leaf element
                leaf->items = allocate_memory(1, int);
                // if memory was allocated successfully
                if (leaf->items)
                {
                        // macro to  initialize a leaf node 'n' from cluster 'c'
                        initialize_leaf(c, leaf, item, length);
                        c->clusters_count++;
                }
                else
                {
                        // print error message to the standard error output and free up space for the tag of the node
                        fprintf(stderr, "Failed to allocate memory for node items.\n");
                        free(leaf->tag);
                        leaf = NULL;
                }
        }
        else
        {
                // print error to stanrad error output
                fprintf(stderr, "Failed to allocate memory for node tag.\n");
                leaf = NULL;
        }
        return leaf;
}

// do not need to use macro any further
#undef initialize_leaf

// Function to add multiple leaf nodes to a cluster
cluster_type *add_leave_nodes(cluster_type *c, labeled_data_type *items)
{
        for (int i = 0; i < c->items_count; ++i)
        {
                if (add_leaf(c, &items[i]))
                        // after each leaf element is added, we need to upgrade neighbours for 
                        // the newly added element
                        update_neighbour_nodes(c, i);
                else
                {
                        // if element wasn't added due to an error, exist the loop and nullify cluster
                        c = NULL;
                        break;
                }
        }
        // return cluster
        return c;
}

// Function to print tags for the data points, which are part of the cluster
void display_cluster_elements(cluster_type *c, int index)
{
        // get node from the cluster 'c' by index
        node_type *curr = &(c->nodes[index]);
        fprintf(stdout, "Items: ");
        // proceed with printing tags only of there are any elements
        if (curr->items_count > 0)
        {
                // display first element without leading space
                fprintf(stdout, "%s", c->nodes[curr->items[0]].tag);
                // and all subsequent ones with the leading space
                for (int i = 1; i < curr->items_count; ++i)
                        fprintf(stdout, ", %s", c->nodes[curr->items[i]].tag);
        }
        // add new line character when all nodes were displayed
        fprintf(stdout, "\n");
}

// Function to print full infromation about the node in the cluster 'c'
void display_cluster_node(cluster_type *c, int index)
{
        // get node from the cluster by index
        node_type *curr = &(c->nodes[index]);
        // Print height and centroid of the node
        fprintf(stdout, "Node[%d]: height = %d, centroid = (%6.2f, %6.2f)\n",
                index, curr->height, curr->centroid.x, curr->centroid.y);
        // print tag for data nodes and 'merged' otheriwse
        if (curr->tag)
                fprintf(stdout, "\t Leaf: %s\n\t", curr->tag);
        else
                fprintf(stdout, "\t Merged: %d, %d\n\t",
                        curr->merged[0], curr->merged[1]);
        // display tags for all child elements
        display_cluster_elements(c, index);
        // display node's neighbours
        fprintf(stdout, "\t Neighbours: ");
        neighbour_type *t = curr->neighbours;
        // while there is next neighbour exists
        while (t)
        {
                fprintf(stdout, "\n\t\t%2d: %5.3f", t->point, t->dist);
                t = t->next;
        }
        fprintf(stdout, "\n");
}

// Function to merge node 'n_to_merge' into node 'n' inside the cluster 'c'
void merge_items(cluster_type *c, node_type *n, node_type **n_to_merge)
{
        n->type = NODE_MERGED_ROOTS;
        n->is_root = 1;
        n->height = -1;

        int k = 0, idx;
        coord_type centroid = {.x = .0, .y = .0};
        // 
        for (int i = 0; i < 2; ++i)
        {
                node_type *t = n_to_merge[i];
                // reset is root property
                t->is_root = 0; 
                // set height (or level) of the node in the cliter
                if (n->height == -1 || n->height < t->height) {
                        n->height = t->height;
                }
                // refresh indexes of the node elements
                for (int j = 0; j < t->items_count; ++j)
                {
                        idx = t->items[j];
                        n->items[k++] = idx;
                }
                // adjust centroid per node, being merged
                centroid.x += t->items_count * t->centroid.x;
                centroid.y += t->items_count * t->centroid.y;
        }
        // calculate node's centroid after update
        n->centroid.x = centroid.x / k;
        n->centroid.y = centroid.y / k;
        n->height++;
}

// Macro to merge nodes 'to_merge' with node 'n' that has index 'n_index'
#define merge_two_nodes(c, to_merge, n, n_index)                                           \
{       /* add items to the node 'n' */                                                    \
        n->items_count = to_merge[0]->items_count + to_merge[1]->items_count;              \
        n->items = allocate_memory(n->items_count, int);                                   \
        if (n->items)                                                                      \
        {                                                                                  \
                /* Sync nodes. Decrement cluster count and increment nodes count */        \
                merge_items(c, n, to_merge);                                               \
                c->nodes_count++;                                                          \
                c->clusters_count--;                                                       \
                update_neighbour_nodes(c, n_index);                                        \
        }                                                                                  \
        else {                                                                             \
                /* in case when memory in the node 'n' could not be allocated for 'items'*/\
                /* display an error to standard error output */                            \
                fprintf(stderr, "Failed to allocate memory for array of merged items.\n"); \
                free(n->merged);                                                           \
                n = NULL;                                                                  \
        }                                                                                  \
}

// Function to merge two nodes with indexes 'index1' and 'index2' together with their 
// child elements into the cluster 'c'
node_type *merge(cluster_type *c, int index1, int index2)
{
        int index_n = c->nodes_count;
        node_type *n = &(c->nodes[index_n]);
        n->merged = allocate_memory(2, int);
        if (n->merged)
        {
                // memory allocation successfull
                node_type *to_merge[2] = { &(c->nodes[index1]), &(c->nodes[index2]) };
                n->merged[0] = index1;
                n->merged[1] = index2;
                // merge nodes and add them to the cluster 'c' 
                merge_two_nodes(c, to_merge, n, index_n);
        }
        else
        {
                // print error to stanrad error output
                fprintf(stderr, "Failed to allocate memory for array of merged nodes");
                n = NULL;
        }
        return n;
}

// No longer will be using this macro 
#undef merge_two_nodes

// Search for the neighbour, which has a distance less than 'smallest_distance'
void search_neighgbour_by_distance(node_type *nodes, int node_idx, neighbour_type *nn, 
        float *smallest_distance, int *index1, int *index2)
{
        while (nn)
        {
                // if node is root, then check for distance and stop
                if (nodes[nn->point].is_root)
                {
                        // if it is a merge node or neighbour has smaller distance than searched 'smallest_distance'
                        if (*index1 == -1 || nn->dist < *smallest_distance)
                        {
                                // set output ref parameters index1 and index2 
                                *index1 = node_idx;
                                *index2 = nn->point;
                                // set smallest distance to the current one
                                *smallest_distance = nn->dist;
                        }
                        break;
                }
                // if node is not a root element, continue, untill one is found
                nn = nn->next;
        }
}

// Traversing cluster from top and down, finding the two nodes that must be merged
int get_clusters_to_merge(cluster_type *c, int *index1, int *index2)
{
        // initialize shortest distance with 0
        float shortest_distance = 0.0;
        int root_clusters_visited = 0;
        int j = c->nodes_count;
        *index1 = -1;
        // if visited not all nodes, then continue the loop
        while (root_clusters_visited < c->clusters_count)
        {
                node_type *node = &(c->nodes[--j]);

                // process only nodes, that contain data
                if (node->type == NODE_NOT_IN_USE || !node->is_root)
                        continue;
                
                ++root_clusters_visited;
                // on the first pass shortest distance will be initialized by the distance of the first merge node or
                // the smaller distance, if found by traversing the cluster up
                search_neighgbour_by_distance(c->nodes, j, node->neighbours, &shortest_distance, index1, index2);
        }
        return *index1;
}

// Funciton to create all merge nodes in the cluster, organizing its data elements into
// an hierarchical structure
cluster_type *try_merge_clusters(cluster_type *c)
{
        int index1, index2;
        while (c->clusters_count > 1)
        {
                // get indexes of the elements, which must be merged
                if (get_clusters_to_merge(c, &index1, &index2) != -1) {
                        // perform a merge of two nodes
                        merge(c, index1, index2);
                }
        }
        return c;
}

// Macro to create a new cluster 'c' from the leaf nodes 'items' and after adding them,
// organize them in an hierarchical structure
#define initiazlie_cluster(c, items_count, items)                 \
        {       /* create distance matrix from all items, */      \
                /* which we will be iused to create a cluster */  \
                c->dist = create_dist_matrix(items_count, items); \
                /* if distance matrix was not created, stop */    \
                /* the execution and free up used memory */       \
                if (!c->dist)                                     \
                        goto free_memory;                         \
                c->items_count = items_count;                     \
                c->nodes_count = 0;                               \
                c->clusters_count = 0;                            \
                /* if adding leaf nodes failed, then stop */      \
                /* the execution and free up used memory */       \
                if (add_leave_nodes(c, items))                    \
                        try_merge_clusters(c);                    \
                else                                              \
                        goto free_memory;                         \
        }                                                         \


// An entry point, used to create a cluster based on the data elements, and organize them
// in a hierarchical structure in cluster
cluster_type *agglomerate(int items_count, labeled_data_type *items)
{
        cluster_type *c = allocate_memory(1, cluster_type);
        if (c)
        {
                // memory allocation for cluster was successfull
                c->nodes = allocate_memory(2 * items_count - 1, node_type);
                if (c->nodes) {
                        // memory allocation for nodes was successfull - continue with creating cluster
                        initiazlie_cluster(c, items_count, items);
                }
                else
                {
                        // if memory allocation for nodes failed, then print error to stanrad error output
                        fprintf(stderr, "Failed to allocate memory for cluster nodes");
                        goto free_memory;
                }
        }
        else
                // if memory allocation failed, then print error to stanrad error output
                fprintf(stderr, "Failed to allocate memory for cluster");
        goto ready;

free_memory:
        // creation of cluster failed - do memory cleanup and return a NULL value
        free_cluster_memory(c);
        c = NULL;

ready:
        // return newly created cluster
        return c;
}

// Macro will not be used beyond this point
#undef initiazlie_cluster

// Display all child nodes for the the merged node with index 'i', ignoring 
// first 'n_ignore' elements (data nodes).
// Function returns integer - count of nodes, which were displayed.
int display_children(cluster_type *c, int i, int n_ignore)
{
        node_type *node = &(c->nodes[i]);
        int root_count = 0;
        // display only merged nodes
        if (node->type == NODE_MERGED_ROOTS)
        {
                // show both children of the merged node
                for (int j = 0; j < 2; j++)
                {
                        int t = node->merged[j];
                        // unless reached index element with index 'n_ignore'
                        if (t < n_ignore)
                        {
                                display_cluster_elements(c, t);
                                ++root_count;
                        }
                }
        }
        return root_count;
}

// Display 'k' nodes from the cluster c
void display_clusters(cluster_type *c, int k)
{
        // if requested to display 0 element, do not do anything
        if (k < 1)
                return;
        // if requested to display more elements that cluster has, set 
        // k to the number of elements in the cluster
        if (k > c->items_count) 
        {
                k = c->items_count;
        }

        // use variable i as an index of the currently iterated element
        int i = c->nodes_count - 1;
        int root_count = 0;
        // calculate count of nodes that need to be skipped
        int nodes_to_ignore = c->nodes_count - k + 1;
        while (k)
        {
                if (i < nodes_to_ignore)
                {
                        // show parent merged nodes
                        display_cluster_elements(c, i);
                        root_count = 1;
                }
                else 
                {
                        // show nodes from the cluster 'c'  from node with index 'x', 
                        // starting from the element, following node with index 'nodes_to_ignore'
                        root_count = display_children(c, i, nodes_to_ignore);
                }
                // decrement counter by the number of displayed nodes
                k -= root_count;
                --i;
        }
}

// Display all emenents of the cluser 'c'
void display_cluster(cluster_type *c)
{
        for (int i = 0; i < c->nodes_count; ++i)
        {
                // Display full information about each node
                display_cluster_node(c, i);
        }
}

// Load 'len' elements of labeled data from file 'file' into memory, adding them to the array 'items'
int load_data(int len, labeled_data_type *items, FILE *file)
{
        // 'len' number times read data element from file
        for (int i = 0; i < len; ++i)
        {
                labeled_data_type *t = &(items[i]);
                // if line does matches required format, then proceed with reading the next one
                if (fscanf(file, "%[^|]| %10f %10f\n", t->tag, &(t->coordinate.x), &(t->coordinate.y)))
                {
                        continue;
                }
                // if line format did not match, then print error to stanrad error output and stop loading
                fprintf(stderr, "Failed to read item line from file.\n");
                // return number of elements, read from file in case of an error
                return i;
        }
        // if no errors, results will match to the number of requested lines to be loaded
        return len;
}

// Read all data elements from file into array of labeled data 'ld'
int load_data_from_file(labeled_data_type **ld, FILE *file)
{
        int len, r;
        // first line in the file will be a number of data elements to be loaded
        r = fscanf(file, "%5d\n", &len);
        // if number of elements is 0, then stop loading
        if (r == 0)
        {
                // print error to stanrad error output, indicating to number of lines could not be loaded 
                // from the file
                fprintf(stderr, "Failed to read number of lines from file (first line in the file \
                        must be a positive integer, indicating number of data elements in the file).\n");
                return 0;
        }
        if (len)
        {
                // prepare for loading data from file by requsting memory for 'len' number of labeled data nodes
                *ld = allocate_memory(len, labeled_data_type);
                if (*ld)
                {
                        // if memory was successfully allocated, then load it into 'ld' array
                        if (load_data(len, *ld, file) != len) 
                        {
                                // if error happened while loading data, free up used memory
                                free(ld);
                        }
                }
                else
                {
                        // print error to stanrad error output, indicating the root cause of the error 
                        fprintf(stderr, "Failed to allocate memory for labeled data elements array");
                }
        }
        return len;
}

// set delegate function for calculating distances between nodes, based on specified value 'selected_type'
//  in the command line argument. Single link function is used by default.
void specify_link(char selected_type)
{
        switch (selected_type)
        {
                case LINKAGE_AVG:
                        dist_float_ptr = avg_link;
                        break;
                case LINKAGE_COMPLETE:
                        dist_float_ptr = complete_link;
                        break;
                case LINKAGE_CENTROID:
                        dist_float_ptr = centr_link;
                        break;
                case LINKAGE_SINGLE:
                default:
                        dist_float_ptr = single_link;
        }
}

// Read number of elements and then all labeled data elements from file, stored on hard drive
// in the file with name, specified in the 'filename' parameter
int read_file(labeled_data_type **ld, const char *filename)
{
        int size = 0;
        // open file in read mode
        FILE *f = fopen(filename, "r");
        if (f)
        {
                // load all data from file
                size = load_data_from_file(ld, f);
                // close file handler
                fclose(f);
        }
        else
        {
                // print error to stanrad error output, indicating that file could not be opened
                fprintf(stderr, "Failed to open input file %s.\n", filename);
        }
        // number of elements, loaded from file
        return size;
}

// Entry point to the application
int main(int argc, char **argv)
{
        // expected number of arguments is always 4
        if (argc != 4)
        {
                // if actual number of arguments not maches the expected couunt,
                // print usage instructions to the stanrad error output
                fprintf(stderr, "Usage:\n\t\t %s <file name> <nummber of clusters> <link type>\n", argv[0]);
                // and stop execution
                exit(1);
        }
        else
        {
                // create a pointer to the labeled data array
                labeled_data_type *ld = NULL;
                // read data from file into the array
                int items_count = read_file(&ld, argv[1]);
                // set delegate, which will be used for calculating distance between nodes, based on the 
                // selected option, passed as a command line argument
                specify_link(argv[3][0]);

                // if there were any elements, loaded from file, start processing
                if (items_count)
                {
                        // create a cluster, based on labeled data elements 'ld'
                        cluster_type *c = agglomerate(items_count, ld);
                        // dispose of memory, taken up by data elements, as their values were already copied into 
                        // the cluster nodes
                        free(ld);

                        // if cluster was created successfully
                        if (c)
                        {
                                // for debugging purposes output of the cluster information can be skipped by not 
                                // defining PRINT_HIERARCHY
                                #ifdef PRINT_HIERARCHY

                                // display hierarchy of elements in the cluster
                                fprintf(stdout, "Hierarchy:\n\n----------------------------------------\n\n");
                                display_cluster(c);
                                
                                #endif

                                // for debugging purposes output of all tags in each cluster can be skipped by not 
                                // defining PRINT_HIERARCHY
                                #ifdef PRINT_CLUSTERS

                                // display labels of data elements inside of each of the clusters
                                int k = atoi(argv[2]);
                                fprintf(stdout, "\n\n%d Clusters:\n\n----------------------------------------\n\n", k);
                                display_clusters(c, k);

                                #endif

                                // free up memory, used by the cluster and all of it's nodes
                                free_cluster_memory(c);
                        }
                }
        }
        // return code of successful execution
        return 0;
}