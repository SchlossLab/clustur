# Create Sparse Matrix

Given a list of i indexes, j indexes, and distances values, we can
create a sparse distance matrix for you. Each vector must have the same
size.

## Usage

``` r
create_sparse_matrix(i_index, j_index, distances)
```

## Arguments

- i_index:

  A list of i indexes, must be numeric

- j_index:

  A list of j indexes, must be numeric

- distances:

  A list of the distance at the i and j index

## Value

a \`dgTMatrix\` from the \`Matrix\` library.

## Examples

``` r
 i_values <- as.integer(1:100)
 j_values <- as.integer(sample(1:100, 100, TRUE))
 x_values <- as.numeric(runif(100, 0, 1))
 s_matrix <- create_sparse_matrix(i_values, j_values, x_values)
```
