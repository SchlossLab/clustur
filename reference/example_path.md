# Example Path

This function was created as a helper function to generate file paths to
our internal data. You should use this function if you want to follow
along with the example, or interact with the data

## Usage

``` r
example_path(file = NULL)
```

## Arguments

- file:

  The file name of the data; leave as NULL (default) to get full list of
  example files

## Value

the path to the file as a \`character\` or a vector of \`character\`
giving example filenames if \`fill = NULL\`.

## Examples

``` r
example_path("amazon_phylip.dist")
#> [1] "/home/runner/work/_temp/Library/clustur/extdata/amazon_phylip.dist"
example_path()
#> [1] "amazon.full.count_table"   "amazon.sparse.count_table"
#> [3] "amazon_column.dist"        "amazon_phylip.dist"       
```
