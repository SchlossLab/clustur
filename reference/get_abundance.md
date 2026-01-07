# Get Shared

GetShared returns the generated abundance \`data.frame\` from the
\`cluster()\` function

## Usage

``` r
get_abundance(cluster_data)
```

## Arguments

- cluster_data:

  The output from the \`cluster()\` function.

## Value

a shared data.frame

## Examples

``` r
cutoff <- 0.2
count_table <- read_count(example_path("amazon.full.count_table"))
distance_data <- read_dist(example_path("amazon_column.dist"),
                           count_table, cutoff, FALSE)
df_clusters <- cluster(distance_data, cutoff, method = "opticlust")
shared <- get_abundance(df_clusters)
```
