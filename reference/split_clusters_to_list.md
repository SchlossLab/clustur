# Split Clusters to List

\`split_clusters_to_list()\` will extract clusters from the cluster
generated \`data.frame\`. It will then turn those clusters into a list.
This allows users to more easily visualize their data.

## Usage

``` r
split_clusters_to_list(cluster)
```

## Arguments

- cluster:

  The output generated from the \`cluster()\` function.

## Value

a named \`list\` of clusters.

## Examples

``` r
 cutoff <- 0.2
 count_table <- read_count(example_path("amazon.full.count_table"))
 distance_data <- read_dist(example_path("amazon_column.dist"),
                            count_table, cutoff, FALSE)
 cluster_results <- cluster(distance_data, cutoff, method = "opticlust")

 cluster_list <- split_clusters_to_list(cluster_results)
```
