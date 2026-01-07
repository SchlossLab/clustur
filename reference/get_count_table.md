# Get Count Table

This function returns the count table that was used to generate the
distance object.

## Usage

``` r
get_count_table(distance_object)
```

## Arguments

- distance_object:

  The output from the \`read.dist()\` function.

## Value

a count_table \`data.frame\`.

## Examples

``` r
cutoff <- 0.2
count_table <- read_count(example_path("amazon.full.count_table"))
distance_data <- read_dist(example_path("amazon_column.dist"),
                           count_table, cutoff, FALSE)
count_table <- get_count_table(distance_data)
```
