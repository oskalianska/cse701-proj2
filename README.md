# C++ project - Hierarchical clustering 

1. [Summary](#1-summary)
2. [Implementation details](#2-implementation-details)
3. [Showcase](#3-showcase)
4. [References](#references)

## 1. Summary

This project aims to implement a hierarchical clustering algorithm. Hierarchical cluster analysis is an algorithm that groups similar items into clusters (groups). 

The final result is a set of clusters, where each cluster differs from each other cluster, and the items inside each cluster are similar.

The problem's main idea is running an algorithm with different linkage criteria (single, complete, average, and centroid linkage). Also, the number of clusters that we would like to have at the end of running it.

## 2. Implementation details

There are two types of hierarchical clustering: divisive (top-down) and agglomerative (bottom-up). In this program, I have used the second type - agglomerative clustering. It starts with assumption that each data point is a separate cluster. Then the two closest clusters are joined into one cluster. The next nearest clusters are grouped, and this process continues until there is only one cluster containing the entire data set.

So, my main algorithm was:

1. Compute the proximity matrix 
2. Let each data point be a cluster 
3. Repeat 
	4. Merge the two closest clusters 
	5. Update the proximity matrix 
6. Until only a single cluster remains

Key operation is the computation of the nearness of two clusters.

### Measures of distance (similarity)

The distance between two clusters has been computed based on the length of the straight line drawn from one cluster to another. This is commonly referred to as the Euclidean distance. The most common definition is with Euclidean distance, minimizing the Sum of Squared Error (SSE) function.

### Linkage Criteria

After selecting a distance metric, it is necessary to determine from what point we compute the distance. 

In single-link, we merge in each step the two clusters whose two closest members have the smallest distance. It can sometimes produce clusters where observations in different clusters are closer together than to observations within their own clusters. These clusters can appear spread-out.

In complete-link (or complete linkage) hierarchical clustering, we choose the cluster pair whose merge has the smallest diameter. This method usually produces tighter clusters than single-linkage, but these tight clusters can end up very close together. Along with average-linkage, it is one of the more popular distance metrics.

Mean, or average-linkage is computed between the center of the clusters. Average-linkage and complete-linkage are the two most popular distance metrics in hierarchical clustering.

The centroid link is the distance between two clusters is defined as the distance between the centroid for group 1 and the centroid for group 2. As the centroids move with new observations, the smaller clusters may be more similar to the new, more massive cluster than to their individual clusters, causing an inversion in the dendrogram. This problem doesn't arise in the other linkage methods because the clusters being merged will always be more similar to themselves than to the new larger cluster.

## 3. Showcase

To run the program folloing command line arguments must be specified:

1. File name with data of you are working
2. How many clusters you want to have in results
3. What type of linkage criteria you want to use ("a" - average, "t" - centroid, "c" - complete, "s" - single)

For example, if one wants to use a hierarchical clustering on sample_data.txt input file, creating 2 clusters with single-linkage, the following command is required:

`$ ./clustering sample_data.txt 2 s`

The choice of linkage method entirely depends on you and there is no hard and fast method that will always give you good results. Different linkage methods lead to different clusters.

Example of my data you can find in attachment to the letter `sample_data.txt`. It's a simple txt file with 10,000 rows. First line indicates number of elements to load from the file. For the test run I have used 3,000 elements, specifiying it in the very first line of the file. However, it can be changed to load more elements, up to 10,000. Every row has a tag along with two parametrs: x and y. 

For instance:

```
EL_1| 12 74.8
EL_2| 5 18.8
EL_3| 9 2.4
EL_4| 27 71.2
EL_5| 60 41.9
EL_6| 24 88.7
EL_7| 91 93.7
```

This data was generated randomly, but in real life it has many cases of uses. 

An example of out put result of the program will looks like:

```
CLUSTER HIERARCHY
--------------------
Node 0 - height: 0, centroid: (12.000, 74.800)
	Leaf: EL_1
	Items: EL_1
	Neighbours: 
Node 1 - height: 0, centroid: (5.000, 18.800)
	Leaf: EL_2
	Items: EL_2
	Neighbours: 
		 0: 56.436
```

As you can see, hierarchical clustering helps us make predictions and show the relationship between sets and groups. In general, this program is an example of a tool to do hierarchical clustering analysis with different parameters.

## References

* https://towardsdatascience.com/introduction-hierarchical-clustering-d3066c6b560e
* https://nlp.stanford.edu/IR-book/html/htmledition/single-link-and-complete-link-clustering-1.html
* https://www.datanovia.com/en/lessons/agglomerative-hierarchical-clustering/
* https://www.displayr.com/what-is-hierarchical-clustering/#:~:text=Hierarchical%20clustering%2C%20also%20known%20as,broadly%20similar%20to%20each%20other.
