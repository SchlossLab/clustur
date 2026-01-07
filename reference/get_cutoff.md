# Get Cutoff

Returns the distance cutoff of the cluster object from the \`cluster()\`
function

## Usage

``` r
get_cutoff(cluster_data)
```

## Arguments

- cluster_data:

  The output from the \`cluster()\` function.

## Value

the cutoff value as a \`dbl\`

## Examples

``` r
cutoff <- 0.2
count_table <- read_count(example_path("amazon.full.count_table"))
distance_data <- read_dist(example_path("amazon_column.dist"),
                           count_table, cutoff, FALSE)
df_clusters <- cluster(distance_data, cutoff, method = "opticlust")
cutoff <- get_cutoff(df_clusters)
```
