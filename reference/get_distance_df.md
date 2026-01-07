# Get Distance Data Frame

This function will generate a \`data.frame\` that contains the distances
of all the indexes.

## Usage

``` r
get_distance_df(distance_object)
```

## Arguments

- distance_object:

  The output from the \`read.dist()\` function.

## Value

a distance \`data.frame\`.

## Examples

``` r
cutoff <- 0.2
count_table <- read_count(example_path("amazon.full.count_table"))
distance_data <- read_dist(example_path("amazon_column.dist"),
                                        count_table, cutoff, FALSE)
count_table <- get_count_table(distance_data)
```
