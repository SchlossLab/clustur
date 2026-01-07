# Get Clusters

GetClusters returns a \`data.frame\` of the generated clusters from the
\`cluster()\` function.

## Usage

``` r
get_bins(cluster_data)
```

## Arguments

- cluster_data:

  The output from the \`cluster()\` function.

## Value

the created cluster \`data.frame\`.

## Examples

``` r
cutoff <- 0.2
count_table <- read_count(example_path("amazon.full.count_table"))
distance_data <- read_dist(example_path("amazon_column.dist"),
                           count_table, cutoff, FALSE)
df_clusters <- cluster(distance_data, cutoff, method = "opticlust")
clusters <- get_bins(df_clusters)
```
