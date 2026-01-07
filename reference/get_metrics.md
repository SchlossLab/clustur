# Get Metrics

GetMetrics returns the generated metrics \`data.frame\` from the
\`cluster()\` function.

## Usage

``` r
get_metrics(cluster_data)
```

## Arguments

- cluster_data:

  The output from the \`cluster()\` function.

## Value

a list of metric data.frames

## Examples

``` r
cutoff <- 0.2
count_table <- read_count(example_path("amazon.full.count_table"))
distance_data <- read_dist(example_path("amazon_column.dist"),
                           count_table, cutoff, FALSE)
df_clusters <- cluster(distance_data, cutoff, method = "opticlust")
list_of_metrics <- get_metrics(df_clusters)
```
